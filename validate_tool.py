#!/usr/bin/env python3
"""Measure what the shipped tool does, by planting swaps we know the answer to.

A QC tool that cannot state its own false-alarm rate is not usable, so this
plants known swaps in both validation cohorts, runs the same code path the CLI
runs, and reports how many it caught and how many clean runs it disturbed.

    python validate_tool.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from proteotrace.core import normalise, pair_table, similarity, threshold_from_pairs
from proteotrace.report import run_report

COHORTS = {
    "PXD001064": {
        "matrix": "data/PXD001064_protein_matrix.tsv",
        "samples": "examples/PXD001064_samples.tsv",
        "note": "twin cohort, 2 visits per donor: the hardest case",
    },
    "PXD004242": {
        "matrix": "data/PXD004242/PXD004242_protein_matrix.tsv",
        "samples": "examples/PXD004242_samples.tsv",
        "note": "7 visits per donor, no twins: the typical case",
    },
}


def plant_swaps(samples: pd.DataFrame, n_pairs: int, rng) -> tuple[pd.DataFrame, list[str]]:
    """Swap donor labels between runs on the same plate, as a real mix-up would."""
    table, swapped = samples.copy(), []
    for _, block in table.groupby("batch"):
        if len(swapped) >= n_pairs * 2:
            break
        idx = list(block.index)
        rng.shuffle(idx)
        for a, b in zip(idx[::2], idx[1::2]):
            if len(swapped) >= n_pairs * 2:
                break
            if table.loc[a, "donor"] != table.loc[b, "donor"]:
                table.loc[a, "donor"], table.loc[b, "donor"] = (
                    table.loc[b, "donor"], table.loc[a, "donor"])
                swapped += [table.loc[a, "run"], table.loc[b, "run"]]
    return table, swapped


def evaluate(matrix: pd.DataFrame, samples: pd.DataFrame, margin: float,
             specificity: float, n_pairs: int, repeats: int) -> dict:
    sensitivity, false_alarm = [], []
    for repeat in range(repeats):
        planted, swapped = plant_swaps(samples, n_pairs, np.random.default_rng(100 + repeat))
        annotation = planted.set_index("run")
        normalised = normalise(matrix, annotation["batch"])
        scores = similarity(normalised)
        pairs = pair_table(scores, annotation["donor"], annotation["batch"])
        match, rival = threshold_from_pairs(pairs, specificity)
        report = run_report(scores, annotation["donor"], annotation["batch"],
                            match, margin, rival_threshold=rival)
        report["planted"] = report["run"].isin(swapped)
        flagged = report["verdict"] == "swap_suspected"
        assessable = report["verdict"] != "no_reference"
        sensitivity.append((flagged & report["planted"]).sum()
                           / max((report["planted"] & assessable).sum(), 1))
        false_alarm.append((flagged & ~report["planted"]).sum()
                           / max((~report["planted"] & assessable).sum(), 1))
    return {"sensitivity": float(np.mean(sensitivity)),
            "false_alarm_rate": float(np.mean(false_alarm)),
            "repeats": repeats, "planted_pairs_per_repeat": n_pairs}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--margin", type=float, default=0.20)
    parser.add_argument("--specificity", type=float, default=0.99)
    parser.add_argument("--pairs", type=int, default=12)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--out", type=Path, default=Path("data/tool_validation.json"))
    args = parser.parse_args()

    results = {"margin": args.margin, "specificity": args.specificity, "cohorts": {}}
    print(f"margin {args.margin}  ·  specificity {args.specificity:.0%}\n")
    print(f"{'cohort':<12}{'runs':>6}{'donors':>8}{'sensitivity':>13}{'false alarms':>14}")
    for name, spec in COHORTS.items():
        matrix = pd.read_csv(spec["matrix"], sep="\t", index_col=0)
        samples = pd.read_csv(spec["samples"], sep="\t")
        shared = matrix.index.intersection(samples["run"])
        matrix = matrix.loc[shared]
        samples = samples[samples["run"].isin(shared)].reset_index(drop=True)
        outcome = evaluate(matrix, samples, args.margin, args.specificity,
                           args.pairs, args.repeats)
        outcome.update({"runs": int(len(samples)),
                        "donors": int(samples["donor"].nunique()),
                        "note": spec["note"]})
        results["cohorts"][name] = outcome
        print(f"{name:<12}{outcome['runs']:>6}{outcome['donors']:>8}"
              f"{outcome['sensitivity']:>12.1%}{outcome['false_alarm_rate']:>14.1%}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
