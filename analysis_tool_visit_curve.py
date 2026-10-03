#!/usr/bin/env python3
"""Swap detection by the shipped tool, against visits per donor (Figure 6B).

Replaces the earlier reference-curve figure, which scored a closed-set,
within-plate-group donor identification rather than the tool (audit F11).
Here each donor in PXD004242 is subsampled to k visits, the subsampled cohort
is run through the tool exactly as `proteotrace check` runs it, and the same
planted-swap accounting as validate_tool.py is applied.

    python analysis_tool_visit_curve.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from proteotrace.core import normalise, similarity
from validate_tool import account, plant_swaps, verdicts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--pairs", type=int, default=12)
    ap.add_argument("--margin", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=20261002)
    ap.add_argument("--out", type=Path, default=Path("data/PXD004242/PXD004242_tool_visit_curve.json"))
    args = ap.parse_args()

    matrix = pd.read_csv("data/PXD004242/PXD004242_protein_matrix.tsv", sep="\t", index_col=0)
    samples = pd.read_csv("examples/PXD004242_samples.tsv", sep="\t", dtype=str)
    samples = samples[samples["run"].isin(matrix.index)].reset_index(drop=True)

    rows = []
    for k in range(2, 8):
        for r in range(args.repeats):
            rng = np.random.default_rng(args.seed + 1000 * k + r)
            keep = []
            for _, block in samples.groupby("donor"):
                if len(block) >= k:
                    keep += list(rng.choice(block.index.to_numpy(), k, replace=False))
            sub = samples.loc[sorted(keep)].reset_index(drop=True)
            annotation = sub.set_index("run")
            scores = similarity(normalise(matrix.loc[sub["run"]], annotation["batch"]))
            planted, swapped = plant_swaps(sub, args.pairs, rng)
            row = account(verdicts(scores, planted, args.margin, 0.99), sub, swapped)
            row.update({"n_visits": k, "repeat": r, "donors": int(sub["donor"].nunique()),
                        "planted_runs": len(swapped)})
            rows.append(row)
        print(f"visits={k} done", flush=True)

    table = pd.DataFrame(rows)
    summary = []
    for k, block in table.groupby("n_visits"):
        entry = {"n_visits": int(k), "donors": int(block["donors"].median()),
                 "planted_runs": int(block["planted_runs"].median())}
        for key in ("sensitivity", "flagged_uninvolved", "flagged_collateral"):
            v = block[key].to_numpy()
            entry[key] = {"mean": float(v.mean()), "q025": float(np.quantile(v, 0.025)),
                          "q975": float(np.quantile(v, 0.975))}
        summary.append(entry)
        print(f"k={k} caught {entry['sensitivity']['mean']:.3f}  uninvolved "
              f"{entry['flagged_uninvolved']['mean']:.3f}  collateral {entry['flagged_collateral']['mean']:.3f}")
    args.out.write_text(json.dumps({"margin": args.margin, "repeats": args.repeats,
                                    "pairs_per_repeat": args.pairs, "rows": summary,
                                    "interpretation": "Means and 2.5-97.5% ranges over "
                                    "repeats of visit subsampling plus swap planting; "
                                    "not confidence intervals."},
                                   indent=2) + "\n", encoding="utf-8")
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
