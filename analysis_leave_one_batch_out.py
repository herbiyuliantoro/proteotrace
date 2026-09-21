#!/usr/bin/env python3
"""Conservative leave-one-acquisition-day-out validation for ProteoTrace.

Each fold estimates feature selection and training standardization only from
non-held-out runs.  The held-out day's median is estimated from its unlabeled
runs because that is the deployment operation of the unsupervised correction;
no donor labels enter any normalization step.  Scores are evaluated only
between held-out runs and reference runs from other days.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from analysis_cross_batch import acquisition, read_metadata


def prepare(matrix_path: Path, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    raw = read_metadata(metadata_path)
    common = [run for run in matrix.index if run in raw]
    matrix = matrix.loc[common]
    metadata = pd.DataFrame([raw[run] for run in common], index=common)
    keep, seen = [], set()
    for run in matrix.index:
        key = (metadata.loc[run, "donor"], metadata.loc[run, "visit"])
        if key in seen:
            continue
        seen.add(key)
        keep.append(run)
    matrix = matrix.loc[keep]
    metadata = metadata.loc[keep]
    parsed = [acquisition(run) for run in matrix.index]
    metadata["batch"] = [x[0] for x in parsed]
    metadata["position"] = [x[1] for x in parsed]
    metadata["day"] = [x[2] for x in parsed]
    return matrix, metadata


def residualize(train: pd.DataFrame, test: pd.DataFrame, train_meta: pd.DataFrame, test_meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    train_positive = train.replace(0, np.nan)
    selected = train_positive.notna().mean(axis=0) >= 0.99
    train_log = np.log2(train_positive.loc[:, selected])
    test_log = np.log2(test.replace(0, np.nan).loc[:, selected])
    train_log = train_log.subtract(train_log.median(axis=1), axis=0)
    test_log = test_log.subtract(test_log.median(axis=1), axis=0)
    train_batch_medians = train_log.groupby(train_meta["batch"]).transform("median")
    test_batch_medians = test_log.groupby(test_meta["batch"]).transform("median")
    train_corr = train_log - train_batch_medians
    test_corr = test_log - test_batch_medians
    mean = train_corr.mean(axis=0)
    std = train_corr.std(axis=0).replace(0, np.nan)
    return (train_corr - mean) / std, (test_corr - mean) / std, int(selected.sum())


def pair_score(left: np.ndarray, right: np.ndarray) -> float:
    mask = np.isfinite(left) & np.isfinite(right)
    if int(mask.sum()) < 30:
        return float("nan")
    x, y = left[mask], right[mask]
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt(np.dot(x, x) * np.dot(y, y))
    return float(np.dot(x, y) / denom) if denom > 0 else float("nan")


def fold_scores(train_res: pd.DataFrame, test_res: pd.DataFrame, train_meta: pd.DataFrame, test_meta: pd.DataFrame, heldout: str) -> pd.DataFrame:
    rows = []
    for test_name in test_res.index:
        left = test_res.loc[test_name].to_numpy(dtype=float)
        for train_name in train_res.index:
            score = pair_score(left, train_res.loc[train_name].to_numpy(dtype=float))
            if not np.isfinite(score):
                continue
            rows.append({
                "heldout_batch": heldout,
                "test_run": test_name,
                "reference_run": train_name,
                "score": score,
                "same_donor": bool(test_meta.loc[test_name, "donor"] == train_meta.loc[train_name, "donor"]),
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=Path("data/PXD001064_protein_matrix.tsv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/PXD001064_sample_info.txt"))
    parser.add_argument("--out-json", type=Path)
    parser.add_argument("--out-tsv", type=Path)
    args = parser.parse_args()
    matrix, metadata = prepare(args.matrix, args.metadata)
    rows = []
    folds = []
    for heldout in sorted(metadata["batch"].unique()):
        test_names = metadata.index[metadata["batch"] == heldout]
        train_names = metadata.index[metadata["batch"] != heldout]
        train_res, test_res, n_features = residualize(
            matrix.loc[train_names], matrix.loc[test_names], metadata.loc[train_names], metadata.loc[test_names]
        )
        fold = fold_scores(train_res, test_res, metadata.loc[train_names], metadata.loc[test_names], heldout)
        folds.append(fold)
        if fold.empty or fold["same_donor"].nunique() < 2:
            auc = None
        else:
            auc = float(roc_auc_score(fold["same_donor"].astype(int), fold["score"]))
        summary = {
            "heldout_batch": heldout,
            "train_runs": int(len(train_names)),
            "test_runs": int(len(test_names)),
            "selected_features": n_features,
            "pairs": int(len(fold)),
            "same_donor_pairs": int(fold["same_donor"].sum()) if not fold.empty else 0,
            "auc": auc,
        }
        rows.append(summary)
        print(
            f"heldout_batch={heldout} train_runs={summary['train_runs']} test_runs={summary['test_runs']} "
            f"selected_features={n_features} pairs={summary['pairs']} same={summary['same_donor_pairs']} auc={auc if auc is not None else 'NA'}"
        )
    pair_table = pd.concat(folds, ignore_index=True) if folds else pd.DataFrame()
    if not pair_table.empty:
        # The same unordered cross-batch pair is seen from either fold. Average
        # those two predictions before the pooled AUC so it is not counted twice.
        pair_table["pair_key"] = pair_table.apply(lambda r: "|".join(sorted([r["test_run"], r["reference_run"]])), axis=1)
        unique = pair_table.groupby("pair_key", as_index=False).agg(score=("score", "mean"), same_donor=("same_donor", "first"))
        pooled_auc = float(roc_auc_score(unique["same_donor"].astype(int), unique["score"])) if unique["same_donor"].nunique() > 1 else None
    else:
        unique = pair_table
        pooled_auc = None
    fold_aucs = [x["auc"] for x in rows if x["auc"] is not None]
    result = {
        "method": "leave_one_acquisition_day_out",
        "normalization": "training feature selection and standardization; held-out batch median from unlabeled held-out runs",
        "folds": rows,
        "mean_fold_auc": float(np.mean(fold_aucs)) if fold_aucs else None,
        "pooled_unique_cross_batch_pairs": int(len(unique)),
        "pooled_unique_auc": pooled_auc,
    }
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(result, indent=2) + "\n")
    if args.out_tsv:
        args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
        pair_table.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"mean_fold_auc={result['mean_fold_auc']}")
    print(f"pooled_unique_cross_batch_pairs={result['pooled_unique_cross_batch_pairs']} pooled_unique_auc={pooled_auc}")


if __name__ == "__main__":
    main()
