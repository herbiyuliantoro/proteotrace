#!/usr/bin/env python3
"""Export the PXD004242 primary pair scores so figures can show distributions.

The frozen analyses store summary AUCs only, which is enough for a point
estimate but not enough to draw a distribution. This script re-runs the same
prespecified pipeline for one assignment and writes every primary pair with its
raw and corrected score, so a figure can show the overlap the AUC summarises.

The assignment is the deterministic timepoint-modulo-4 rotation, and the pair
set is the primary estimand: different actual plates inside one plate group.
The printed AUCs must match the frozen rotation result (raw 0.8425,
corrected 0.9839); the script fails loudly if they do not.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_pxd004242_assignment_sensitivity import residual_scores, select_rotation
from analysis_pxd004242_layer1 import read_sample_table

EXPECTED = {"raw": 0.8425, "corrected": 0.9839}
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
                        default=Path("data/PXD004242/PXD004242_primary_pairs_rotation.tsv"))
    args = parser.parse_args()

    base = args.root / "data/PXD004242"
    matrix = pd.read_csv(base / "PXD004242_protein_matrix.tsv", sep="\t", index_col=0)
    rows = read_sample_table(base / "PXD004242_biological_samples.tsv", "longitudinal")
    metadata = select_rotation(rows).copy()
    metadata["plate_group"] = metadata["plate"].str[0]
    matrix = matrix.loc[metadata.index]

    raw, corrected, n_features = residual_scores(matrix, metadata)

    n = len(metadata)
    ii, jj = np.triu_indices(n, 1)
    donor = metadata["individual"].to_numpy()
    plate = metadata["plate"].to_numpy()
    group = metadata["plate_group"].to_numpy()

    primary = (group[ii] == group[jj]) & (plate[ii] != plate[jj])
    finite = np.isfinite(raw[ii, jj]) & np.isfinite(corrected[ii, jj])
    keep = primary & finite
    ii, jj = ii[keep], jj[keep]

    pairs = pd.DataFrame({
        "left": metadata.index.to_numpy()[ii],
        "right": metadata.index.to_numpy()[jj],
        "plate_group": group[ii],
        "same_donor": donor[ii] == donor[jj],
        "score_raw": raw[ii, jj],
        "score_corrected": corrected[ii, jj],
    })

    observed = {"raw": auc(pairs["same_donor"].to_numpy(), pairs["score_raw"].to_numpy()),
                "corrected": auc(pairs["same_donor"].to_numpy(), pairs["score_corrected"].to_numpy())}
    for key, value in observed.items():
        print(f"{key:>10s} AUC {value:.4f}  (frozen {EXPECTED[key]:.4f})")
    drift = {k: abs(observed[k] - EXPECTED[k]) for k in EXPECTED}
    if max(drift.values()) > TOLERANCE:
        raise SystemExit(f"AUC does not reproduce the frozen rotation result: {drift}")

    out = args.root / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(out, sep="\t", index=False)
    print(f"features={n_features} runs={n} pairs={len(pairs)} "
          f"same_donor={int(pairs['same_donor'].sum())}")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
