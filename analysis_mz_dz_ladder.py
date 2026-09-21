#!/usr/bin/env python3
"""MZ/DZ/unrelated similarity ladder for the PXD001064 pilot."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from analysis_cross_batch import acquisition, pair_table, read_metadata


def prepare(matrix_path: Path, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    raw = read_metadata(metadata_path)
    common = [run for run in matrix.index if run in raw]
    matrix = matrix.loc[common]
    metadata = pd.DataFrame([raw[run] for run in common], index=common)

    keep, seen = [], set()
    for run in matrix.index:
        key = (metadata.loc[run, "donor"], metadata.loc[run, "visit"])
        if key in seen:
            continue
        seen.add(key)
        keep.append(run)
    matrix = matrix.loc[keep]
    metadata = metadata.loc[keep]

    parsed = [acquisition(run) for run in matrix.index]
    metadata["batch"] = [x[0] for x in parsed]
    metadata["position"] = [x[1] for x in parsed]
    metadata["day"] = [x[2] for x in parsed]

    positive = matrix.replace(0, np.nan)
    selected = positive.notna().mean(axis=0) >= 0.99
    log_matrix = np.log2(positive.loc[:, selected])
    log_matrix = log_matrix.subtract(log_matrix.median(axis=1), axis=0)
    batch_medians = log_matrix.groupby(metadata["batch"]).transform("median")
    residual = log_matrix - batch_medians
    residual = (residual - residual.mean(axis=0)) / residual.std(axis=0).replace(0, np.nan)
    return residual.T.corr(min_periods=30), metadata


def relation(left: pd.Series, right: pd.Series) -> str:
    if left["donor"] == right["donor"]:
        return "same_donor"
    left_pair = re.match(r"(Twin pair-\d+)-[AB]$", str(left["donor"]))
    right_pair = re.match(r"(Twin pair-\d+)-[AB]$", str(right["donor"]))
    if left_pair and right_pair and left_pair.group(1) == right_pair.group(1):
        return "MZ_twins" if left["zygosity"] == "MZ" else "DZ_twins"
    return "unrelated"


def summarize(table: pd.DataFrame, label: str) -> None:
    values = table["score"]
    print(
        f"{label}: n={len(values)} median={values.median():.4f} "
        f"q05={values.quantile(.05):.4f} q25={values.quantile(.25):.4f} "
        f"q75={values.quantile(.75):.4f} q95={values.quantile(.95):.4f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=Path("data/PXD001064_protein_matrix.tsv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/PXD001064_sample_info.txt"))
    args = parser.parse_args()

    scores, metadata = prepare(args.matrix, args.metadata)
    pairs = pair_table(scores, metadata)
    pairs["relation"] = [relation(metadata.loc[a], metadata.loc[b]) for a, b in zip(pairs["left"], pairs["right"])]
    pairs = pairs.loc[~pairs["same_batch"]].copy()

    print(f"cross_batch_pairs={len(pairs)}")
    print("relation_counts", pairs["relation"].value_counts().to_dict())
    for label in ["same_donor", "MZ_twins", "DZ_twins", "unrelated"]:
        subset = pairs.loc[pairs["relation"] == label]
        summarize(subset, label)

    order = [("same_donor", "MZ_twins"), ("MZ_twins", "DZ_twins"), ("DZ_twins", "unrelated"), ("same_donor", "unrelated")]
    for left, right in order:
        x = pairs.loc[pairs["relation"] == left, "score"]
        y = pairs.loc[pairs["relation"] == right, "score"]
        u, p_value = mannwhitneyu(x, y, alternative="two-sided")
        probability = u / (len(x) * len(y))
        print(f"compare_{left}_vs_{right}: P(left>right)={probability:.4f} p={p_value:.3g}")


if __name__ == "__main__":
    main()
