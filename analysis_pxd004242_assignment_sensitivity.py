#!/usr/bin/env python3
"""Assignment-aware layer-1 sensitivity for PXD004242.

Each draw keeps exactly one run per biological sample. Primary pairs are
restricted to different actual plates inside the same plate group. This
removes the donor-group/plate-group nesting from the negative class while
retaining cross-plate same-donor positives.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from analysis_pxd004242_layer1 import TIME_ORDER, read_sample_table, select_runs


SPECIFICITIES = (0.95, 0.99)


def corr_scores(values: pd.DataFrame) -> np.ndarray:
    # pandas' pairwise correlation preserves the existing min-periods contract.
    return values.T.corr(min_periods=30).to_numpy(dtype=float)


def residual_scores(matrix: pd.DataFrame, metadata: pd.DataFrame, feature_columns: list[str] | None = None):
    positive = matrix.replace(0, np.nan)
    if feature_columns is None:
        selected = positive.notna().mean(axis=0) >= 0.99
        log_values = np.log2(positive.loc[:, selected])
        n_features = int(selected.sum())
    else:
        log_values = np.log2(positive.loc[:, feature_columns])
        n_features = len(feature_columns)
    log_values = log_values.subtract(log_values.median(axis=1), axis=0)
    raw_corr = corr_scores(log_values)

    batch_medians = log_values.groupby(metadata["plate"]).transform("median")
    residual = log_values - batch_medians
    residual = (residual - residual.mean(axis=0)) / residual.std(axis=0).replace(0, np.nan)
    corrected_corr = corr_scores(residual)
    return raw_corr, corrected_corr, n_features


def metrics_for_scores(scores: np.ndarray, metadata: pd.DataFrame) -> dict[str, object]:
    n = len(metadata)
    ii, jj = np.triu_indices(n, 1)
    values = scores[ii, jj]
    finite = np.isfinite(values)
    values = values[finite]
    ii, jj = ii[finite], jj[finite]
    donor = metadata["individual"].to_numpy()
    plate = metadata["plate"].to_numpy()
    group = metadata["plate_group"].to_numpy()
    same_donor = donor[ii] == donor[jj]
    same_plate = plate[ii] == plate[jj]
    same_group = group[ii] == group[jj]
    primary = same_group & ~same_plate
    same_group_all = same_group
    all_cross_plate = ~same_plate

    def auc(mask):
        y = same_donor[mask]
        return float(roc_auc_score(y.astype(int), values[mask])) if y.size and np.unique(y).size == 2 else None

    result = {
        "all_pairs": {"n": int(values.size), "same": int(same_donor.sum()), "auc": auc(np.ones(values.size, dtype=bool))},
        "within_group_all_pairs": {"n": int(same_group_all.sum()), "same": int(same_donor[same_group_all].sum()), "auc": auc(same_group_all)},
        "cross_plate_all_groups": {"n": int(all_cross_plate.sum()), "same": int(same_donor[all_cross_plate].sum()), "auc": auc(all_cross_plate)},
        "primary_within_group_cross_plate": {
            "n": int(primary.sum()),
            "same": int(same_donor[primary].sum()),
            "different": int((~same_donor[primary]).sum()),
            "auc": auc(primary),
        },
    }
    primary_scores = values[primary]
    primary_same = same_donor[primary]
    for specificity in SPECIFICITIES:
        negatives = primary_scores[~primary_same]
        positives = primary_scores[primary_same]
        key = f"specificity_{int(specificity * 100)}"
        if not len(negatives) or not len(positives):
            result["primary_within_group_cross_plate"][key] = None
        else:
            threshold = float(np.quantile(negatives, specificity))
            result["primary_within_group_cross_plate"][key] = {
                "threshold": threshold,
                "sensitivity": float(np.mean(positives >= threshold)),
                "n_positive": int(len(positives)),
                "n_negative": int(len(negatives)),
            }
    return result


def select_rotation(rows):
    by_sample = {}
    for row in rows:
        by_sample.setdefault(row["sample_id"], []).append(row)
    selected = []
    for sample_id in sorted(by_sample):
        candidates = by_sample[sample_id]
        target = ((TIME_ORDER[candidates[0]["sampling_time"]] - 1) % len(candidates)) + 1
        matches = [r for r in candidates if int(r["technical_replicate"]) == target]
        if len(matches) != 1:
            raise ValueError(f"rotation failed for {sample_id}")
        selected.append(matches[0])
    return pd.DataFrame(selected).set_index("run")


def one_assignment(matrix: pd.DataFrame, metadata: pd.DataFrame, feature_columns: list[str] | None = None) -> dict[str, object]:
    metadata = metadata.copy()
    metadata["plate_group"] = metadata["plate"].str[0]
    raw, corrected, n_features = residual_scores(matrix.loc[metadata.index], metadata, feature_columns)
    result = {
        "n_runs": int(len(metadata)),
        "n_donors": int(metadata["individual"].nunique()),
        "n_actual_plates": int(metadata["plate"].nunique()),
        "n_plate_groups": int(metadata["plate_group"].nunique()),
        "selected_features": n_features,
        "raw": metrics_for_scores(raw, metadata),
        "corrected": metrics_for_scores(corrected, metadata),
    }
    return result


def summarize_draws(draws: list[dict[str, object]]) -> dict[str, object]:
    def get(path, draw):
        x = draw
        for key in path:
            x = x[key]
        return x
    paths = {
        "corrected_primary_auc": ["corrected", "primary_within_group_cross_plate", "auc"],
        "raw_primary_auc": ["raw", "primary_within_group_cross_plate", "auc"],
        "corrected_primary_same": ["corrected", "primary_within_group_cross_plate", "same"],
        "corrected_primary_n": ["corrected", "primary_within_group_cross_plate", "n"],
        "corrected_sens_95": ["corrected", "primary_within_group_cross_plate", "specificity_95", "sensitivity"],
        "corrected_sens_99": ["corrected", "primary_within_group_cross_plate", "specificity_99", "sensitivity"],
    }
    summary = {"n_draws": len(draws)}
    for name, path in paths.items():
        vals = np.asarray([get(path, d) for d in draws if get(path, d) is not None], dtype=float)
        if vals.size:
            summary[name] = {
                "median": float(np.median(vals)),
                "q025": float(np.quantile(vals, 0.025)),
                "q975": float(np.quantile(vals, 0.975)),
                "min": float(np.min(vals)),
                "max": float(np.max(vals)),
            }
        else:
            summary[name] = None
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--sample-table", type=Path, required=True)
    ap.add_argument("--mode", choices=["rotation", "random"], required=True)
    ap.add_argument("--n-draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260919)
    ap.add_argument("--feature-mode", choices=["adaptive_99", "top273_global", "random273_global"], default="adaptive_99")
    ap.add_argument("--feature-seed", type=int, default=20260919)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    args = ap.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    rows = read_sample_table(args.sample_table, "longitudinal")
    feature_columns = None
    if args.feature_mode != "adaptive_99":
        available_runs = [r["run"] for r in rows if r["run"] in matrix.index]
        coverage = matrix.loc[available_runs].replace(0, np.nan).notna().mean(axis=0)
        ranked = coverage.sort_values(ascending=False, kind="mergesort").index.tolist()
        if args.feature_mode == "top273_global":
            feature_columns = ranked[:273]
        else:
            rng_features = np.random.default_rng(args.feature_seed)
            feature_columns = sorted(rng_features.choice(matrix.columns.to_numpy(), size=273, replace=False).tolist())
    draws = []
    if args.mode == "rotation":
        metadata = select_rotation(rows)
        draws.append({"draw": 0, "selection": "timepoint_mod_4", **one_assignment(matrix, metadata, feature_columns)})
    else:
        for i in range(args.n_draws):
            draw_seed = args.seed + i
            metadata = select_runs(rows, "random", None, draw_seed)
            draws.append({"draw": i, "seed": draw_seed, **one_assignment(matrix, metadata, feature_columns)})
            if (i + 1) % 20 == 0 or i == 0:
                print(f"completed_draws={i + 1}/{args.n_draws}")

    summary = {
        "method": "one_run_assignment_sensitivity",
        "mode": args.mode,
        "feature_mode": args.feature_mode,
        "feature_seed": args.feature_seed if args.feature_mode == "random273_global" else None,
        "feature_count": len(feature_columns) if feature_columns is not None else None,
        "negative_restriction": "different actual plates within the same plate group (plate-group is the leading digit)",
        "cohort": "longitudinal_only",
        "draw_summary": summarize_draws(draws),
        "draws": draws,
        "interpretation": "Quantiles describe run-selection variability, not a statistical confidence interval. Plate-group restriction removes the donor-nested-group negative confound.",
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rows_tsv = []
    for d in draws:
        p = d["corrected"]["primary_within_group_cross_plate"]
        r = d["raw"]["primary_within_group_cross_plate"]
        rows_tsv.append({
            "draw": d["draw"], "seed": d.get("seed", ""), "n_runs": d["n_runs"],
            "n_actual_plates": d["n_actual_plates"], "selected_features": d["selected_features"],
            "raw_auc": r["auc"], "corrected_auc": p["auc"], "n_pairs": p["n"],
            "same_donor_pairs": p["same"], "different_donor_pairs": p["different"],
            "corrected_sens_95": p["specificity_95"]["sensitivity"] if p["specificity_95"] else "",
            "corrected_sens_99": p["specificity_99"]["sensitivity"] if p["specificity_99"] else "",
        })
    pd.DataFrame(rows_tsv).to_csv(args.out_tsv, sep="\t", index=False)
    print(json.dumps(summary["draw_summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
