#!/usr/bin/env python3
"""Cross-batch confounding check for the PXD001064 ProteoTrace pilot.

The matrix is an OpenSWATH target-rank-1 protein-group table produced by the
pilot extraction. No donor labels are used to select features or normalize
the scores. Acquisition batch is inferred from the run filename date and is
used only for an unsupervised per-batch median correction.
"""
from __future__ import annotations

import argparse
import csv
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score


META_KEY = "Twin Pair Information (A and B represent two twins in one pair, see TableS1 in the publication)"


def read_metadata(path: Path) -> dict[str, dict[str, str | int]]:
    text = path.read_text(errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    result: dict[str, dict[str, str | int]] = {}
    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        if row.get("Visit") not in {"1", "2"}:
            continue
        result[row["Raw file names"]] = {
            "donor": row[META_KEY],
            "visit": int(row["Visit"]),
            "zygosity": row["Zygosity"],
        }
    return result


def acquisition(run: str) -> tuple[str, int, int]:
    match = re.search(r"L(\d{6})_(\d+)", run)
    if not match:
        raise ValueError(f"cannot parse acquisition date/order from {run}")
    date_text, number_text = match.groups()
    day = datetime.strptime(date_text, "%y%m%d")
    number = int(number_text)
    # The filename gives a day and within-day run number, not a clock time.
    position = day.toordinal() * 1000 + number
    return date_text, position, day.toordinal()


def pair_table(scores: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    records = []
    names = list(scores.index)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            score = scores.loc[left, right]
            if not np.isfinite(score):
                continue
            records.append(
                {
                    "left": left,
                    "right": right,
                    "score": float(score),
                    "same_donor": metadata.loc[left, "donor"] == metadata.loc[right, "donor"],
                    "same_batch": metadata.loc[left, "batch"] == metadata.loc[right, "batch"],
                    "date_gap_days": abs(int(metadata.loc[left, "day"]) - int(metadata.loc[right, "day"])),
                    "run_gap": abs(int(metadata.loc[left, "position"]) - int(metadata.loc[right, "position"])),
                }
            )
    return pd.DataFrame(records)


def summarize(table: pd.DataFrame, label: str) -> None:
    if table.empty or table["same_donor"].nunique() < 2:
        print(f"{label}: unavailable (n={len(table)})")
        return
    auc = roc_auc_score(table["same_donor"].astype(int), table["score"])
    same = table.loc[table["same_donor"], "score"]
    different = table.loc[~table["same_donor"], "score"]
    print(
        f"{label}: n={len(table)} same={len(same)} different={len(different)} "
        f"auc={auc:.4f} same_median={same.median():.4f} different_median={different.median():.4f}"
    )
    if len(same) and len(different):
        threshold = different.quantile(0.95)
        print(
            f"{label}: threshold_95pct_specificity={threshold:.4f} "
            f"same_pair_sensitivity={np.mean(same >= threshold):.4f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=Path("/tmp/PXD001064_protein_matrix.tsv"))
    parser.add_argument("--metadata", type=Path, default=Path("/tmp/PXD001064_sample_info.txt"))
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    metadata_raw = read_metadata(args.metadata)
    common = [run for run in matrix.index if run in metadata_raw]
    matrix = matrix.loc[common]
    metadata = pd.DataFrame([metadata_raw[run] for run in common], index=common)

    # The deposited metadata has two duplicated donor/visit entries. Keep one
    # run per donor and visit so the denominator matches the published design.
    keep = []
    seen = set()
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
    # Baseline without batch correction.
    raw_scores = log_matrix.T.corr(min_periods=30)
    raw_pairs = pair_table(raw_scores, metadata)
    # Unsupervised batch correction: per-protein median within acquisition day.
    batch_medians = log_matrix.groupby(metadata["batch"]).transform("median")
    residual = log_matrix - batch_medians
    residual = (residual - residual.mean(axis=0)) / residual.std(axis=0).replace(0, np.nan)
    scores = residual.T.corr(min_periods=30)
    pairs = pair_table(scores, metadata)

    print(f"runs_after_dedup={len(matrix)} donors={metadata['donor'].nunique()} complete_pairs={sum(metadata['visit'].eq(1) & metadata['donor'].isin(metadata.loc[metadata['visit'].eq(2), 'donor']))}")
    print(f"protein_groups_selected={int(selected.sum())} acquisition_batches={metadata['batch'].nunique()}")
    summarize(raw_pairs.loc[~raw_pairs["same_batch"]], "raw_cross_batch_pairs")
    summarize(pairs, "all_pairs")
    summarize(pairs.loc[~pairs["same_batch"]], "cross_batch_pairs")
    summarize(pairs.loc[pairs["date_gap_days"] > 0], "cross_day_pairs")

    for subset_name, subset in [("same_donor", pairs.loc[pairs["same_donor"]]), ("different_donor", pairs.loc[~pairs["same_donor"]])]:
        if len(subset) < 3:
            continue
        rho, p_value = spearmanr(subset["score"], subset["date_gap_days"])
        print(f"score_vs_date_gap_{subset_name}: rho={rho:.4f} p={p_value:.3g} n={len(subset)}")
    u, p_value = mannwhitneyu(
        pairs.loc[pairs["same_donor"], "date_gap_days"],
        pairs.loc[~pairs["same_donor"], "date_gap_days"],
        alternative="two-sided",
    )
    print(
        f"date_gap_same_median={pairs.loc[pairs['same_donor'], 'date_gap_days'].median():.1f} "
        f"different_median={pairs.loc[~pairs['same_donor'], 'date_gap_days'].median():.1f} "
        f"mannwhitney_p={p_value:.3g}"
    )


if __name__ == "__main__":
    main()
