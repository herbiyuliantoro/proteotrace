#!/usr/bin/env python3
"""Leave-one-plate-out batch-transfer check for PXD004242 layer 1."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from analysis_pxd004242_layer1 import read_sample_table, select_runs


def pair_score(left: np.ndarray, right: np.ndarray) -> float:
    mask = np.isfinite(left) & np.isfinite(right)
    if int(mask.sum()) < 30:
        return float("nan")
    x, y = left[mask], right[mask]
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt(np.dot(x, x) * np.dot(y, y))
    return float(np.dot(x, y) / denom) if denom > 0 else float("nan")


def residualize(train: pd.DataFrame, test: pd.DataFrame, train_meta: pd.DataFrame, test_meta: pd.DataFrame):
    train_positive = train.replace(0, np.nan)
    selected = train_positive.notna().mean(axis=0) >= 0.99
    train_log = np.log2(train_positive.loc[:, selected])
    test_log = np.log2(test.replace(0, np.nan).loc[:, selected])
    train_log = train_log.subtract(train_log.median(axis=1), axis=0)
    test_log = test_log.subtract(test_log.median(axis=1), axis=0)
    train_batch_medians = train_log.groupby(train_meta["plate"]).transform("median")
    test_batch_medians = test_log.groupby(test_meta["plate"]).transform("median")
    train_corr = train_log - train_batch_medians
    test_corr = test_log - test_batch_medians
    mean = train_corr.mean(axis=0)
    std = train_corr.std(axis=0).replace(0, np.nan)
    return (train_corr - mean) / std, (test_corr - mean) / std, int(selected.sum())


def fold_scores(train_res, test_res, train_meta, test_meta, heldout):
    rows = []
    for test_name in test_res.index:
        left = test_res.loc[test_name].to_numpy(dtype=float)
        for train_name in train_res.index:
            score = pair_score(left, train_res.loc[train_name].to_numpy(dtype=float))
            if not np.isfinite(score):
                continue
            rows.append({
                "heldout_plate": heldout,
                "test_run": test_name,
                "reference_run": train_name,
                "score": score,
                "same_donor": bool(test_meta.loc[test_name, "individual"] == train_meta.loc[train_name, "individual"]),
            })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--sample-table", type=Path, required=True)
    ap.add_argument("--selection", choices=["fixed", "rotate_by_time", "random"], default="fixed")
    ap.add_argument("--technical-replicate", type=int, choices=[1, 2, 3, 4])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-pairs", type=Path)
    args = ap.parse_args()
    if args.selection == "fixed" and args.technical_replicate is None:
        ap.error("--technical-replicate is required for --selection fixed")

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    rows = read_sample_table(args.sample_table, "longitudinal")
    metadata = select_runs(rows, args.selection, args.technical_replicate, args.seed)
    matrix = matrix.loc[metadata.index]
    fold_summaries = []
    folds = []
    for heldout in sorted(metadata["plate"].unique()):
        test_names = metadata.index[metadata["plate"] == heldout]
        train_names = metadata.index[metadata["plate"] != heldout]
        train_res, test_res, n_features = residualize(matrix.loc[train_names], matrix.loc[test_names], metadata.loc[train_names], metadata.loc[test_names])
        fold = fold_scores(train_res, test_res, metadata.loc[train_names], metadata.loc[test_names], heldout)
        folds.append(fold)
        auc = None if fold.empty or fold["same_donor"].nunique() < 2 else float(roc_auc_score(fold["same_donor"].astype(int), fold["score"]))
        fold_summaries.append({
            "heldout_plate": heldout,
            "train_runs": int(len(train_names)),
            "test_runs": int(len(test_names)),
            "selected_features": n_features,
            "pairs": int(len(fold)),
            "same_donor_pairs": int(fold["same_donor"].sum()) if not fold.empty else 0,
            "auc": auc,
        })
        print(f"heldout_plate={heldout} train_runs={len(train_names)} test_runs={len(test_names)} selected_features={n_features} pairs={len(fold)} same={int(fold['same_donor'].sum()) if not fold.empty else 0} auc={auc if auc is not None else 'NA'}")
    pair_table = pd.concat(folds, ignore_index=True) if folds else pd.DataFrame()
    if not pair_table.empty:
        pair_table["pair_key"] = pair_table.apply(lambda r: "|".join(sorted([r["test_run"], r["reference_run"]])), axis=1)
        unique = pair_table.groupby("pair_key", as_index=False).agg(score=("score", "mean"), same_donor=("same_donor", "first"))
        pooled_auc = float(roc_auc_score(unique["same_donor"].astype(int), unique["score"])) if unique["same_donor"].nunique() > 1 else None
    else:
        unique = pair_table
        pooled_auc = None
    fold_aucs = [x["auc"] for x in fold_summaries if x["auc"] is not None]
    result = {
        "method": "leave_one_plate_out",
        "normalization": "training feature selection and standardization; held-out plate median from unlabeled held-out runs",
        "selection": args.selection,
        "technical_replicate": args.technical_replicate,
        "seed": args.seed if args.selection == "random" else None,
        "cohort": "longitudinal_only",
        "n_runs": int(len(metadata)),
        "n_donors": int(metadata["individual"].nunique()),
        "n_plates": int(metadata["plate"].nunique()),
        "folds": fold_summaries,
        "mean_fold_auc": float(np.mean(fold_aucs)) if fold_aucs else None,
        "pooled_unique_cross_plate_pairs": int(len(unique)),
        "pooled_unique_auc": pooled_auc,
        "interpretation": "This is plate-correction transfer, not donor-held-out validation; donors remain assigned to plate families.",
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.out_pairs:
        args.out_pairs.parent.mkdir(parents=True, exist_ok=True)
        pair_table.to_csv(args.out_pairs, sep="\t", index=False)
    print(f"mean_fold_auc={result['mean_fold_auc']}")
    print(f"pooled_unique_cross_plate_pairs={result['pooled_unique_cross_plate_pairs']} pooled_unique_auc={result['pooled_unique_auc']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
