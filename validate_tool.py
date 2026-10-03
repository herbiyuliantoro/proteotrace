#!/usr/bin/env python3
"""Measure what the shipped tool does, by planting swaps we know the answer to.

A QC tool that cannot state its own false-alarm rate is not usable, so this
plants known swaps in both validation cohorts, runs the same code path the CLI
runs, and reports how many it caught and how many runs it disturbed.

Design, revised after the pre-preprint audit (docs/AUDIT_REPORT.md, F2):

* Swaps are reciprocal label exchanges between two runs on the same batch,
  which is exactly the data a tube swap after plating produces. Pairs are
  drawn uniformly over *all* batches (the first version filled its quota from
  the alphabetically first batch), and no donor is touched by two swaps.
* "Flagged" runs are split three ways. A swap between donors A and B also
  leaves A's and B's *other* runs compared against a wrong reference; with two
  visits per donor those are flagged by construction. They are reported as
  collateral, separately from runs whose donor was not involved at all.
* The clean cohort (nothing planted) is reported too, including `review`.

The similarity matrix does not depend on donor labels, so each cohort is
scored once and only the labels change between repeats.

    python validate_tool.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from proteotrace.core import normalise, pair_table, similarity, threshold_from_pairs
from proteotrace.report import VERDICTS, run_report

COHORTS = {
    "PXD001064": {
        "matrix": "data/PXD001064_protein_matrix.tsv",
        "samples": "examples/PXD001064_samples.tsv",
        "note": "twin cohort, 2 visits per donor: the hardest case",
    },
    "PXD004242": {
        "matrix": "data/PXD004242/PXD004242_protein_matrix.tsv",
        "samples": "examples/PXD004242_samples.tsv",
        "note": "6-7 visits per donor, no twins: the typical case",
    },
}
MARGINS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30)


def plant_swaps(samples: pd.DataFrame, n_pairs: int, rng) -> tuple[pd.DataFrame, list[str]]:
    """Exchange donor labels between same-batch runs, pairs drawn over all batches."""
    table = samples.copy()
    candidates = []
    for _, block in table.groupby("batch"):
        idx = block.index.to_numpy()
        donors = block["donor"].to_numpy()
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                if donors[a] != donors[b]:
                    candidates.append((idx[a], idx[b]))
    touched, swapped = set(), []
    for k in rng.permutation(len(candidates)):
        a, b = candidates[k]
        da, db = table.at[a, "donor"], table.at[b, "donor"]
        if da in touched or db in touched:
            continue
        table.at[a, "donor"], table.at[b, "donor"] = db, da
        touched |= {da, db}
        swapped += [table.at[a, "run"], table.at[b, "run"]]
        if len(swapped) >= 2 * n_pairs:
            break
    return table, swapped


def verdicts(scores, samples, margin, specificity):
    annotation = samples.set_index("run")
    pairs = pair_table(scores, annotation["donor"], annotation["batch"])
    match, rival = threshold_from_pairs(pairs, specificity)
    report = run_report(scores, annotation["donor"], annotation["batch"],
                        match, margin, rival_threshold=rival)
    return report.set_index("run").loc[samples["run"]]


def account(report, true_samples, swapped) -> dict[str, float]:
    true_donor = true_samples.set_index("run").loc[report.index, "donor"]
    involved = set(true_samples.set_index("run").loc[swapped, "donor"])
    planted = report.index.isin(swapped)
    collateral = true_donor.isin(involved).to_numpy() & ~planted
    uninvolved = ~planted & ~collateral
    flagged = (report["verdict"] == "swap_suspected").to_numpy()
    assessable = report["verdict"].isin(["ok", "review", "swap_suspected"]).to_numpy()

    def rate(mask):
        return float((flagged & mask & assessable).sum() / max((mask & assessable).sum(), 1))

    return {"sensitivity": rate(planted), "flagged_uninvolved": rate(uninvolved),
            "flagged_collateral": rate(collateral), "flagged_any_unplanted": rate(~planted)}


def evaluate(matrix, samples, specificity, n_pairs, repeats, seed) -> dict:
    annotation = samples.set_index("run")
    scores = similarity(normalise(matrix, annotation["batch"]))
    clean = verdicts(scores, samples, 0.20, specificity)["verdict"].value_counts()
    out = {"clean_cohort_default_margin": {v: int(clean.get(v, 0)) for v in VERDICTS},
           "margins": {}}
    rows = []
    for r in range(repeats):
        planted, swapped = plant_swaps(samples, n_pairs, np.random.default_rng(seed + r))
        for margin in MARGINS:
            row = account(verdicts(scores, planted, margin, specificity), samples, swapped)
            row["margin"] = margin
            rows.append(row)
    table = pd.DataFrame(rows)
    for margin, block in table.groupby("margin"):
        out["margins"][f"{margin:.2f}"] = {
            key: {"mean": float(block[key].mean()),
                  "mc_se": float(block[key].std(ddof=1) / np.sqrt(len(block)))}
            for key in ("sensitivity", "flagged_uninvolved", "flagged_collateral",
                        "flagged_any_unplanted")}
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--specificity", type=float, default=0.99)
    parser.add_argument("--pairs", type=int, default=12)
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--seed", type=int, default=100)
    parser.add_argument("--out", type=Path, default=Path("data/tool_validation.json"))
    args = parser.parse_args()

    results = {"specificity": args.specificity, "planted_pairs_per_repeat": args.pairs,
               "repeats": args.repeats, "seed": args.seed,
               "design": "reciprocal same-batch swaps drawn uniformly over all batches; "
                         "no donor in two swaps", "cohorts": {}}
    for name, spec in COHORTS.items():
        matrix = pd.read_csv(spec["matrix"], sep="\t", index_col=0)
        samples = pd.read_csv(spec["samples"], sep="\t", dtype=str)
        samples = samples[samples["run"].isin(matrix.index)].reset_index(drop=True)
        matrix = matrix.loc[samples["run"]]
        outcome = evaluate(matrix, samples, args.specificity, args.pairs, args.repeats, args.seed)
        outcome.update({"runs": int(len(samples)), "donors": int(samples["donor"].nunique()),
                        "note": spec["note"]})
        results["cohorts"][name] = outcome

        clean = outcome["clean_cohort_default_margin"]
        print(f"\n{name}: {outcome['runs']} runs, {outcome['donors']} donors  ({spec['note']})")
        print("  clean cohort at margin 0.20: " + ", ".join(f"{k} {v}" for k, v in clean.items()))
        print(f"  {'margin':>6} {'caught':>8} {'uninvolved':>11} {'collateral':>11}")
        for margin, m in outcome["margins"].items():
            print(f"  {margin:>6} {m['sensitivity']['mean']:>8.1%} "
                  f"{m['flagged_uninvolved']['mean']:>11.1%} {m['flagged_collateral']['mean']:>11.1%}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
