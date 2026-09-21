#!/usr/bin/env python3
"""Simulate metadata swaps under realistic plate constraints.

A swap is represented as two selected runs whose assigned donor labels are
exchanged. Detection is evaluated by leave-two-out donor matching: the query is
matched to donor groups using similarity scores, and a swap is detected when the
best-scoring donor differs from the assigned donor label. The two swapped runs
are excluded from all reference groups to avoid a direct self/counterpart leak.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_pxd004242_layer1 import read_sample_table, select_runs
from analysis_pxd004242_assignment_sensitivity import select_rotation


def corrected_scores(matrix: pd.DataFrame, metadata: pd.DataFrame) -> tuple[np.ndarray, int]:
    positive = matrix.replace(0, np.nan)
    selected = positive.notna().mean(axis=0) >= 0.99
    values = np.log2(positive.loc[:, selected])
    values = values.subtract(values.median(axis=1), axis=0)
    batch = values.groupby(metadata["plate"]).transform("median")
    values = values - batch
    values = (values - values.mean(axis=0)) / values.std(axis=0).replace(0, np.nan)
    scores = values.T.corr(min_periods=30).to_numpy(dtype=float)
    return scores, int(selected.sum())


def donor_scores(scores, metadata, query_idx: int, excluded: set[int], candidate_group: str | None):
    donor = metadata["individual"].to_numpy()
    groups = metadata["plate_group"].to_numpy()
    candidates = sorted(set(donor[i] for i in range(len(donor)) if candidate_group is None or groups[i] == candidate_group))
    out = {}
    for d in candidates:
        members = [i for i in range(len(donor)) if donor[i] == d and i not in excluded and i != query_idx]
        vals = [scores[query_idx, i] for i in members if np.isfinite(scores[query_idx, i])]
        out[d] = float(np.median(vals)) if vals else float("nan")
    return out


def best_donor(score_map):
    finite = [(v, d) for d, v in score_map.items() if np.isfinite(v)]
    if not finite:
        return None, float("nan")
    value, donor = max(finite, key=lambda x: (x[0], x[1]))
    return donor, float(value)


def make_pairs(metadata: pd.DataFrame, constraint: str, rng: np.random.Generator, n_swaps: int):
    plate = metadata["plate"].to_numpy()
    donor = metadata["individual"].to_numpy()
    group = metadata["plate_group"].to_numpy()
    pairs = []
    for i in range(len(metadata)):
        for j in range(i + 1, len(metadata)):
            if donor[i] == donor[j]:
                continue
            if constraint == "within_plate" and plate[i] != plate[j]:
                continue
            if constraint == "within_group" and group[i] != group[j]:
                continue
            pairs.append((i, j))
    if not pairs:
        raise ValueError(f"no swap pairs for constraint={constraint}")
    take = min(n_swaps, len(pairs))
    chosen = rng.choice(len(pairs), size=take, replace=False)
    return [pairs[int(k)] for k in chosen]


def run_simulation(matrix, metadata, constraint, n_swaps, seed):
    scores, n_features = corrected_scores(matrix.loc[metadata.index], metadata)
    metadata = metadata.copy()
    metadata["plate_group"] = metadata["plate"].str[0]
    donor = metadata["individual"].to_numpy()
    group = metadata["plate_group"].to_numpy()
    rng = np.random.default_rng(seed)
    pairs = make_pairs(metadata, constraint, rng, n_swaps)
    rows = []
    for swap_id, (i, j) in enumerate(pairs):
        excluded = {i, j}
        for q, assigned in ((i, donor[j]), (j, donor[i])):
            candidate_group = group[q] if constraint != "global" else None
            baseline_map = donor_scores(scores, metadata, q, excluded, candidate_group)
            top_before, score_before = best_donor(baseline_map)
            assigned_score = baseline_map.get(assigned, float("nan"))
            top_after = top_before
            detected = top_after is not None and top_after != assigned
            rows.append({
                "swap_id": swap_id,
                "query_run": metadata.index[q],
                "plate": metadata.iloc[q]["plate"],
                "plate_group": group[q],
                "true_donor": donor[q],
                "assigned_donor_after_swap": assigned,
                "top_donor_from_profile": top_after,
                "baseline_top_donor": top_before,
                "baseline_correct": bool(top_before == donor[q]),
                "swap_detected": bool(detected),
                "true_donor_score": float(baseline_map.get(donor[q], np.nan)),
                "assigned_donor_score": float(assigned_score),
                "top_score": float(score_before),
                "margin_true_minus_assigned": float(baseline_map.get(donor[q], np.nan) - assigned_score) if np.isfinite(assigned_score) else float("nan"),
            })
    table = pd.DataFrame(rows)
    summary = {
        "constraint": constraint,
        "seed": seed,
        "n_swaps": len(pairs),
        "n_queries": len(table),
        "selected_features": n_features,
        "candidate_group_policy": "same leading-digit plate group" if constraint != "global" else "all groups",
        "baseline_top1_accuracy": float(table["baseline_correct"].mean()),
        "swap_detection_sensitivity": float(table["swap_detected"].mean()),
        "swap_detection_when_baseline_true": float(table.loc[table["baseline_correct"], "swap_detected"].mean()),
        "median_margin_true_minus_assigned": float(table["margin_true_minus_assigned"].median()),
        "plate_counts": table["plate"].value_counts().sort_index().to_dict(),
    }
    return summary, table


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--sample-table", type=Path, required=True)
    ap.add_argument("--selection", choices=["rotation", "random"], default="random")
    ap.add_argument("--seed", type=int, default=20260919)
    ap.add_argument("--n-swaps", type=int, default=1000)
    ap.add_argument("--constraint", choices=["within_plate", "within_group", "global"], default="within_plate")
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    args = ap.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    rows = read_sample_table(args.sample_table, "longitudinal")
    metadata = select_rotation(rows) if args.selection == "rotation" else select_runs(rows, "random", None, args.seed)
    metadata["plate_group"] = metadata["plate"].str[0]
    summary, table = run_simulation(matrix, metadata, args.constraint, args.n_swaps, args.seed)
    summary.update({"selection": args.selection, "cohort": "longitudinal_only", "interpretation": "Detection is top donor disagreement after exchanging two labels; swapped query/reference runs are left out of donor medians."})
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    table.to_csv(args.out_tsv, sep="\t", index=False)
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
