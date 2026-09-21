#!/usr/bin/env python3
"""Estimate within-plate swap detection as a function of donor references.

The original within-plate stress test uses every available longitudinal run in
each donor's reference profile.  That makes its detection rate depend on the
number of visits per donor. This script keeps the same selected runs and the
same 1,000 within-plate swaps, but caps candidate donor profiles at a specified
number of total visits. The query run and the two swapped runs are always
excluded from the reference profiles.

The curve is indexed by total visits per donor. A query donor therefore has at
most ``n_visits - 1`` references, while another candidate donor has at most
``n_visits`` references. The resulting quantiles describe reference-subset
resampling, not a statistical confidence interval and not a validation estimate
for an unknown cohort.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_pxd004242_layer1 import read_sample_table, select_runs
from analysis_pxd004242_swap_simulation import corrected_scores, make_pairs


def prepare_query(
    scores: np.ndarray,
    metadata: pd.DataFrame,
    query_idx: int,
    excluded: set[int],
    candidate_group: str | None,
) -> dict[str, object]:
    """Cache donor-specific score vectors for one query.

    The expensive score lookup and donor membership filtering do not depend on
    visit count or the reference-subset repeat, so they are done once per query.
    """
    donor = metadata["individual"].to_numpy()
    groups = metadata["plate_group"].to_numpy()
    candidates = sorted(
        set(
            donor[i]
            for i in range(len(donor))
            if candidate_group is None or groups[i] == candidate_group
        )
    )
    members_by_donor: dict[str, np.ndarray] = {}
    values_by_donor: dict[str, np.ndarray] = {}
    for label in candidates:
        members = np.asarray(
            [
                i
                for i in range(len(donor))
                if donor[i] == label and i not in excluded and i != query_idx
            ],
            dtype=int,
        )
        members_by_donor[label] = members
        values_by_donor[label] = scores[query_idx, members]
    return {
        "labels": candidates,
        "members": members_by_donor,
        "values": values_by_donor,
    }


def donor_profile_score(
    prepared: dict[str, object],
    n_visits: int,
    true_label: str,
    rng: np.random.Generator,
) -> dict[str, float]:
    labels = prepared["labels"]
    members = prepared["members"]
    values_by_donor = prepared["values"]
    counts = {
        label: n_visits - 1 if label == true_label else n_visits
        for label in labels
    }
    eligible = [label for label in labels if len(members[label]) >= counts[label]]
    if not eligible:
        return {}
    max_len = max(len(values_by_donor[label]) for label in eligible)
    values = np.full((len(eligible), max_len), np.nan, dtype=float)
    for row, label in enumerate(eligible):
        current = np.asarray(values_by_donor[label], dtype=float)
        values[row, : len(current)] = current
    # Random keys select the required number of distinct valid references per
    # donor. Invalid padded cells are assigned an infinite key.
    keys = rng.random(values.shape)
    keys[~np.isfinite(values)] = np.inf
    selected_values = np.full((len(eligible), n_visits), np.nan, dtype=float)
    for row, label in enumerate(eligible):
        count = counts[label]
        selected = np.argpartition(keys[row], kth=count - 1)[:count]
        selected_values[row, :count] = values[row, selected]
    medians = np.nanmedian(selected_values, axis=1)
    return {label: float(value) for label, value in zip(eligible, medians) if np.isfinite(value)}


def best_donor(score_map: dict[str, float]) -> str | None:
    finite = [(value, label) for label, value in score_map.items() if np.isfinite(value)]
    if not finite:
        return None
    return max(finite, key=lambda item: (item[0], item[1]))[1]


def evaluate_curve(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    n_swaps: int,
    n_repeats: int,
    seed: int,
    max_k: int,
    constraint: str,
) -> tuple[dict[str, object], pd.DataFrame]:
    metadata = metadata.copy()
    metadata["plate_group"] = metadata["plate"].str[0]
    scores, n_features = corrected_scores(matrix.loc[metadata.index], metadata)
    rng_pairs = np.random.default_rng(seed)
    pairs = make_pairs(metadata, constraint, rng_pairs, n_swaps)
    donor = metadata["individual"].to_numpy()
    group = metadata["plate_group"].to_numpy()

    query_rows: list[dict[str, object]] = []
    for i, j in pairs:
        excluded = {i, j}
        for query_idx, assigned in ((i, donor[j]), (j, donor[i])):
            query_rows.append(
                {
                    "assigned": assigned,
                    "true": donor[query_idx],
                    "prepared": prepare_query(
                        scores,
                        metadata,
                        query_idx,
                        excluded,
                        group[query_idx],
                    ),
                }
            )

    rows: list[dict[str, object]] = []
    # A separate deterministic stream per (visit count, repeat) makes the output stable
    # if the number of swaps or repeats is changed in a later rerun.
    for n_visits in range(2, max_k + 1):
        for repeat in range(n_repeats):
            rng = np.random.default_rng(seed + 100_000 * n_visits + repeat)
            baseline = []
            detected = []
            correct_recovery = []
            baseline_and_detected = []
            eligible = 0
            for query in query_rows:
                score_map = donor_profile_score(
                    query["prepared"],
                    n_visits,
                    str(query["true"]),
                    rng,
                )
                if str(query["true"]) not in score_map:
                    # A donor with fewer than n_visits total observations is
                    # not eligible for that point on the curve.
                    continue
                top = best_donor(score_map)
                if top is None:
                    continue
                eligible += 1
                is_baseline_correct = top == query["true"]
                is_detected = top != query["assigned"]
                baseline.append(is_baseline_correct)
                detected.append(is_detected)
                correct_recovery.append(is_baseline_correct and is_detected)
                baseline_and_detected.append(is_detected if is_baseline_correct else np.nan)
            rows.append(
                {
                    "n_total_visits": n_visits,
                    "repeat": repeat,
                    "n_queries": len(query_rows),
                    "n_eligible": eligible,
                    "baseline_top1_accuracy": float(np.mean(baseline)) if baseline else np.nan,
                    "swap_detection_sensitivity": float(np.mean(detected)) if detected else np.nan,
                    "correct_swap_recovery": float(np.mean(correct_recovery)) if correct_recovery else np.nan,
                    "swap_detection_when_baseline_true": float(
                        np.nanmean(baseline_and_detected)
                    )
                    if baseline_and_detected and not np.all(np.isnan(baseline_and_detected))
                    else np.nan,
                }
            )

    table = pd.DataFrame(rows)
    summary_rows = []
    for n_visits, group_df in table.groupby("n_total_visits", sort=True):
        row = {
            "n_total_visits": int(n_visits),
            "n_repeats": int(len(group_df)),
            "n_queries": int(group_df["n_queries"].iloc[0]),
            "n_eligible_median": float(group_df["n_eligible"].median()),
        }
        for column in (
            "baseline_top1_accuracy",
            "swap_detection_sensitivity",
            "correct_swap_recovery",
            "swap_detection_when_baseline_true",
        ):
            values = group_df[column].dropna().to_numpy(dtype=float)
            row[f"{column}_median"] = float(np.median(values))
            row[f"{column}_q025"] = float(np.quantile(values, 0.025))
            row[f"{column}_q975"] = float(np.quantile(values, 0.975))
        summary_rows.append(row)

    summary = {
        "cohort": "longitudinal_only",
        "constraint": constraint,
        "seed": seed,
        "n_swaps": len(pairs),
        "n_queries": len(query_rows),
        "n_repeats": n_repeats,
        "selected_features": n_features,
        "max_total_visits": max_k,
        "candidate_group_policy": "same leading-digit plate group",
        "interpretation": (
            "n_total_visits is the assumed number of visits per donor; the query donor "
            "contributes n_total_visits-1 references and other candidates contribute "
            "n_total_visits references after exclusions. Quantiles describe reference-subset resampling, "
            "not a statistical confidence interval or unknown-cohort validation."
        ),
        "summary": summary_rows,
    }
    return summary, table


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--sample-table", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260919)
    parser.add_argument("--n-swaps", type=int, default=1000)
    parser.add_argument("--n-repeats", type=int, default=50)
    parser.add_argument("--max-visits", type=int, default=7)
    parser.add_argument("--constraint", choices=["within_plate", "within_group"], default="within_plate")
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    rows = read_sample_table(args.sample_table, "longitudinal")
    metadata = select_runs(rows, "random", None, args.seed)
    summary, table = evaluate_curve(
        matrix,
        metadata,
        n_swaps=args.n_swaps,
        n_repeats=args.n_repeats,
        seed=args.seed,
        max_k=args.max_visits,
        constraint=args.constraint,
    )
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    table.to_csv(args.out_tsv, sep="\t", index=False)
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
