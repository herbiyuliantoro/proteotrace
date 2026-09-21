#!/usr/bin/env python3
"""Time-distance analysis of donor similarity in PXD004242."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from analysis_pxd004242_assignment_sensitivity import residual_scores
from analysis_pxd004242_layer1 import read_sample_table, select_runs

TIME_WEEKS = {"-8W": -8, "0W": 0, "4W": 4, "13W": 13, "26W": 26, "39W": 39, "52W": 52}


def one_draw(matrix, metadata):
    metadata = metadata.copy()
    metadata["plate_group"] = metadata["plate"].str[0]
    raw_corr, corrected_corr, n_features = residual_scores(matrix.loc[metadata.index], metadata)
    n = len(metadata)
    ii, jj = np.triu_indices(n, 1)
    values = corrected_corr[ii, jj]
    finite = np.isfinite(values)
    values, ii, jj = values[finite], ii[finite], jj[finite]
    donor = metadata["individual"].to_numpy()
    plate = metadata["plate"].to_numpy()
    group = metadata["plate_group"].to_numpy()
    times = metadata["sampling_time"].to_numpy()
    same = donor[ii] == donor[jj]
    cross_plate = plate[ii] != plate[jj]
    same_group = group[ii] == group[jj]
    primary = same_group & cross_plate
    negatives = values[primary & ~same]
    records = []
    gap_values = np.array([abs(TIME_WEEKS[times[a]] - TIME_WEEKS[times[b]]) for a, b in zip(ii, jj)], dtype=int)
    # All -8W-to-nonnegative pairs are the intervention-crossing category.
    pre_post = np.array([(TIME_WEEKS[times[a]] < 0 and TIME_WEEKS[times[b]] >= 0) or (TIME_WEEKS[times[b]] < 0 and TIME_WEEKS[times[a]] >= 0) for a, b in zip(ii, jj)], dtype=bool)
    post_only = np.array([(TIME_WEEKS[times[a]] >= 0 and TIME_WEEKS[times[b]] >= 0) for a, b in zip(ii, jj)], dtype=bool)
    for label, mask in [("pre_to_post", primary & same & pre_post), ("post_only", primary & same & post_only)]:
        pos = values[mask]
        records.append({
            "category": label,
            "gap": "all",
            "n_same": int(len(pos)),
            "same_median": float(np.median(pos)) if len(pos) else None,
            "auc": float(roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(negatives))], np.r_[pos, negatives])) if len(pos) and len(negatives) else None,
        })
    for gap in sorted(set(gap_values[primary & same])):
        pos = values[primary & same & (gap_values == gap)]
        records.append({
            "category": "gap_weeks",
            "gap": int(gap),
            "n_same": int(len(pos)),
            "same_median": float(np.median(pos)) if len(pos) else None,
            "auc": float(roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(negatives))], np.r_[pos, negatives])) if len(pos) and len(negatives) else None,
        })
    return records, n_features


def summarize(rows):
    frame = pd.DataFrame(rows)
    out = []
    for (category, gap), g in frame.groupby(["category", "gap"], sort=True):
        result = {"category": category, "gap": gap, "n_draws": int(len(g))}
        for col in ["n_same", "same_median", "auc"]:
            vals = pd.to_numeric(g[col], errors="coerce").dropna().to_numpy(dtype=float)
            if len(vals):
                result[col + "_median"] = float(np.median(vals))
                result[col + "_q025"] = float(np.quantile(vals, 0.025))
                result[col + "_q975"] = float(np.quantile(vals, 0.975))
            else:
                result[col + "_median"] = result[col + "_q025"] = result[col + "_q975"] = None
        out.append(result)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--sample-table", type=Path, required=True)
    ap.add_argument("--n-draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260919)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    args = ap.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    rows = read_sample_table(args.sample_table, "longitudinal")
    all_rows = []
    feature_counts = []
    for i in range(args.n_draws):
        seed = args.seed + i
        metadata = select_runs(rows, "random", None, seed)
        records, n_features = one_draw(matrix, metadata)
        for r in records:
            r.update({"draw": i, "seed": seed})
        all_rows.extend(records)
        feature_counts.append(n_features)
        if (i + 1) % 20 == 0 or i == 0:
            print(f"completed_draws={i + 1}/{args.n_draws}")
    summary = {
        "method": "time_distance_same_donor_similarity",
        "n_draws": args.n_draws,
        "cohort": "longitudinal_only",
        "negative_restriction": "different actual plates within the same leading-digit plate group",
        "feature_count_median": float(np.median(feature_counts)),
        "summary": summarize(all_rows),
        "interpretation": "Quantiles reflect random one-run assignment variability. Pre_to_post means a pair crosses -8W to a nonnegative timepoint; post_only uses two nonnegative timepoints.",
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    pd.DataFrame(all_rows).to_csv(args.out_tsv, sep="\t", index=False)
    print(json.dumps(summary["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
