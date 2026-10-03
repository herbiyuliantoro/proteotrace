#!/usr/bin/env python3
"""Export the PXD004242 primary pair scores so figures can show distributions.

The frozen analyses store summary AUCs only, which is enough for a point
estimate but not enough to draw a distribution. This script re-runs the
pipeline for one assignment and writes every primary pair with its raw,
scaling-only and corrected score, so a figure can show the overlap the AUC
summarises.

The assignment is the *median* random draw of analysis_primary_estimand.py:
the draw whose corrected AUC is closest to the median of the 200. (The earlier
timepoint-rotation assignment scored above all 200 random draws, so drawing
the figures from it showed a best case; audit F4.) The pair set is the primary
estimand: same-donor against different-donor pairs on different plates, across
plate groups.

The expected values below were verified with an independent implementation
(audit/audit_lib.py, which reads only the MaxQuant zip and the SDRF). The
script fails loudly if the recomputed AUCs drift from them.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_primary_estimand import log_centred, score
from analysis_pxd004242_layer1 import read_sample_table, select_runs

EXPECTED = {"seed": 20261090, "raw": 0.8389, "scaled": 0.9588, "corrected": 0.9674,
            "features": 150, "pairs": 41912, "same_donor": 666}
TOLERANCE = 5e-4


def auc(labels: np.ndarray, values: np.ndarray) -> float:
    order = np.argsort(values)
    ranked = np.empty(len(order), dtype=float)
    ranked[order] = np.arange(1, len(order) + 1)
    pos, neg = labels, ~labels
    return float((ranked[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * neg.sum()))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--out", type=Path,
                        default=Path("data/PXD004242/PXD004242_primary_pairs_median.tsv"))
    args = parser.parse_args()

    base = args.root / "data/PXD004242"
    summary = json.loads((base / "PXD004242_primary_random200.json").read_text(encoding="utf-8"))
    seed = int(summary["median_draw"]["seed"])
    if seed != EXPECTED["seed"]:
        raise SystemExit(f"median draw is seed {seed}, expected {EXPECTED['seed']}: "
                         "the 200-draw analysis changed; re-verify before exporting")

    matrix = pd.read_csv(base / "PXD004242_protein_matrix.tsv", sep="\t", index_col=0)
    rows = read_sample_table(base / "PXD004242_biological_samples.tsv", "longitudinal")
    metadata = select_runs(rows, "random", None, seed)
    log_values = log_centred(matrix.loc[metadata.index])
    scores = {mode: score(log_values, metadata["plate"], mode)
              for mode in ("raw", "scaled", "corrected")}

    n = len(metadata)
    ii, jj = np.triu_indices(n, 1)
    donor = metadata["individual"].to_numpy()
    plate = metadata["plate"].to_numpy()
    group = np.asarray([p[0] for p in plate])
    keep = plate[ii] != plate[jj]
    for values in scores.values():
        keep &= np.isfinite(values[ii, jj])
    ii, jj = ii[keep], jj[keep]

    pairs = pd.DataFrame({
        "left": metadata.index.to_numpy()[ii],
        "right": metadata.index.to_numpy()[jj],
        "same_plate_group": group[ii] == group[jj],
        "same_donor": donor[ii] == donor[jj],
        **{f"score_{mode}": values[ii, jj] for mode, values in scores.items()},
    })

    labels = pairs["same_donor"].to_numpy()
    drift = {}
    for mode in scores:
        observed = auc(labels, pairs[f"score_{mode}"].to_numpy())
        print(f"{mode:>10s} AUC {observed:.4f}  (expected {EXPECTED[mode]:.4f})")
        drift[mode] = abs(observed - EXPECTED[mode])
    counts = {"features": int(log_values.shape[1]), "pairs": int(len(pairs)),
              "same_donor": int(labels.sum())}
    if max(drift.values()) > TOLERANCE or any(counts[k] != EXPECTED[k] for k in counts):
        raise SystemExit(f"does not reproduce the verified median draw: drift={drift} counts={counts}")

    out = args.root / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(out, sep="\t", index=False)
    print(f"seed={seed} features={counts['features']} runs={n} pairs={counts['pairs']} "
          f"same_donor={counts['same_donor']}")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
