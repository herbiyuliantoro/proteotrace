#!/usr/bin/env python3
"""Loaders shared by the main figures.

Keeping these in one place means every figure reads the same frozen artefacts
and no figure can quietly recompute a number differently from its neighbour.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_cross_batch import acquisition, pair_table, read_metadata
from analysis_mz_dz_ladder import prepare, relation


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# ------------------------------------------------------------- PXD001064
def pxd001064_pairs(root: Path) -> pd.DataFrame:
    """Cross-batch pair table with a relatedness label on every row."""
    scores, metadata = prepare(root / "data/PXD001064_protein_matrix.tsv",
                               root / "data/PXD001064_sample_info.txt")
    pairs = pair_table(scores, metadata)
    pairs["relation"] = [relation(metadata.loc[a], metadata.loc[b])
                         for a, b in zip(pairs["left"], pairs["right"])]
    return pairs.loc[~pairs["same_batch"]].copy()


def pxd001064_raw_vs_corrected(root: Path) -> dict[str, object]:
    """Day counts plus the raw and corrected cross-batch AUC."""
    matrix = pd.read_csv(root / "data/PXD001064_protein_matrix.tsv", sep="\t", index_col=0)
    raw_meta = read_metadata(root / "data/PXD001064_sample_info.txt")
    common = [run for run in matrix.index if run in raw_meta]
    matrix = matrix.loc[common]
    metadata = pd.DataFrame([raw_meta[run] for run in common], index=common)

    keep, seen = [], set()
    for run in matrix.index:
        key = (str(metadata.loc[run, "donor"]), int(metadata.loc[run, "visit"]))
        if key not in seen:
            seen.add(key)
            keep.append(run)
    matrix, metadata = matrix.loc[keep], metadata.loc[keep].copy()
    metadata["batch"] = [acquisition(run)[0] for run in matrix.index]

    positive = matrix.replace(0, np.nan)
    selected = positive.notna().mean(axis=0) >= 0.99
    log_matrix = np.log2(positive.loc[:, selected])
    log_matrix = log_matrix.subtract(log_matrix.median(axis=1), axis=0)
    corrected = log_matrix - log_matrix.groupby(metadata["batch"]).transform("median")
    corrected = (corrected - corrected.mean(axis=0)) / corrected.std(axis=0).replace(0, np.nan)

    def cross_batch(scores: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        labels, values = [], []
        names = list(scores.index)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                if metadata.loc[a, "batch"] == metadata.loc[b, "batch"]:
                    continue
                value = scores.loc[a, b]
                if np.isfinite(value):
                    labels.append(metadata.loc[a, "donor"] == metadata.loc[b, "donor"])
                    values.append(float(value))
        return np.asarray(labels, dtype=bool), np.asarray(values, dtype=float)

    def auc(labels: np.ndarray, values: np.ndarray) -> float:
        order = np.argsort(values)
        ranked = np.empty(len(order), dtype=float)
        ranked[order] = np.arange(1, len(order) + 1)
        pos, neg = labels, ~labels
        return float((ranked[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * neg.sum()))

    raw_labels, raw_values = cross_batch(log_matrix.T.corr(min_periods=30))
    cor_labels, cor_values = cross_batch(corrected.T.corr(min_periods=30))
    return {
        "day_counts": dict(sorted(Counter(metadata["batch"].tolist()).items())),
        "raw": {"labels": raw_labels, "values": raw_values, "auc": auc(raw_labels, raw_values)},
        "corrected": {"labels": cor_labels, "values": cor_values, "auc": auc(cor_labels, cor_values)},
        "n_features": int(selected.sum()),
    }


def pxd001064_illustration(root: Path, n_proteins: int = 14, per_batch: int = 5) -> dict:
    """A small real slice of PXD001064 used to illustrate the computation.

    Returns the log-abundance matrix before and after per-batch median
    subtraction for a handful of proteins and runs drawn from two acquisition
    days, plus one same-donor pair measured on different days.
    """
    matrix = pd.read_csv(root / "data/PXD001064_protein_matrix.tsv", sep="\t", index_col=0)
    raw_meta = read_metadata(root / "data/PXD001064_sample_info.txt")
    common = [run for run in matrix.index if run in raw_meta]
    matrix = matrix.loc[common]
    metadata = pd.DataFrame([raw_meta[run] for run in common], index=common)
    metadata["batch"] = [acquisition(run)[0] for run in matrix.index]

    keep, seen = [], set()
    for run in matrix.index:
        key = (str(metadata.loc[run, "donor"]), int(metadata.loc[run, "visit"]))
        if key not in seen:
            seen.add(key)
            keep.append(run)
    matrix, metadata = matrix.loc[keep], metadata.loc[keep]

    positive = matrix.replace(0, np.nan)
    selected = positive.notna().mean(axis=0) >= 0.99
    log_matrix = np.log2(positive.loc[:, selected])
    log_matrix = log_matrix.subtract(log_matrix.median(axis=1), axis=0)
    corrected = log_matrix - log_matrix.groupby(metadata["batch"]).transform("median")

    # Anchor the example on a *representative* same-donor pair: the cross-batch
    # same-donor pair whose score is closest to the median of that class. The
    # first pair encountered is often an outlier and would mislead.
    standardised = (corrected - corrected.mean(axis=0)) / corrected.std(axis=0).replace(0, np.nan)
    scores = standardised.T.corr(min_periods=30)
    candidates = []
    for donor, block in metadata.groupby("donor"):
        runs_of_donor = list(block.index)
        for i, a in enumerate(runs_of_donor):
            for b in runs_of_donor[i + 1:]:
                if metadata.loc[a, "batch"] == metadata.loc[b, "batch"]:
                    continue
                value = scores.loc[a, b]
                if np.isfinite(value):
                    candidates.append((float(value), a, b))
    if not candidates:
        raise ValueError("no same-donor pair spans two acquisition days")
    target = float(np.median([c[0] for c in candidates]))
    _, left, right = min(candidates, key=lambda c: abs(c[0] - target))
    batch_a = metadata.loc[left, "batch"]
    batch_b = metadata.loc[right, "batch"]

    runs = [left] + [r for r in metadata.index[metadata["batch"] == batch_a] if r != left][:per_batch - 1]
    runs += [right] + [r for r in metadata.index[metadata["batch"] == batch_b] if r != right][:per_batch - 1]
    spread = log_matrix.loc[runs].std(axis=0).sort_values(ascending=False)
    proteins = list(spread.index[:n_proteins])
    # Order rows by their between-batch difference. This is ordering, not
    # selection: the same proteins are shown, arranged so the block structure
    # the correction removes is legible.
    block = log_matrix.loc[runs, proteins]
    first = [r for r in runs if metadata.loc[r, "batch"] == batch_a]
    second = [r for r in runs if metadata.loc[r, "batch"] == batch_b]
    delta = block.loc[first].mean(axis=0) - block.loc[second].mean(axis=0)
    proteins = list(delta.sort_values().index)

    return {
        "raw": log_matrix.loc[runs, proteins],
        "corrected": corrected.loc[runs, proteins],
        "batches": metadata.loc[runs, "batch"].tolist(),
        "pair": (left, right),
        "pair_scores": (standardised.loc[left], standardised.loc[right]),
        "pair_r": float(scores.loc[left, right]),
    }


def operating_points(labels: np.ndarray, values: np.ndarray,
                     specificities=(0.95, 0.99)) -> dict[float, dict[str, float]]:
    negatives, positives = values[~labels], values[labels]
    out = {}
    for sp in specificities:
        threshold = float(np.quantile(negatives, sp))
        out[sp] = {"threshold": threshold, "sensitivity": float(np.mean(positives >= threshold))}
    return out


# ------------------------------------------------------------- PXD004242
def pxd004242_design(root: Path) -> pd.DataFrame:
    """One row per donor: its plate group and the rounds it was run in."""
    samples = pd.read_csv(root / "data/PXD004242/PXD004242_biological_samples.tsv", sep="\t")
    samples = samples.loc[samples["cohort"] == "longitudinal"]
    rows = []
    for donor, block in samples.groupby("individual"):
        plates: set[str] = set()
        for text in block["plates"]:
            plates.update(str(text).split(","))
        groups = sorted({p[0] for p in plates})
        rows.append({"donor": int(donor), "group": groups[0] if len(groups) == 1 else "mixed",
                     "rounds": "".join(sorted(p[-1] for p in plates)),
                     "n_timepoints": int(block.shape[0])})
    return pd.DataFrame(rows).sort_values(["group", "donor"]).reset_index(drop=True)


def pxd004242_pair_scores(root: Path) -> dict[str, dict[str, object]]:
    """Primary PXD004242 pairs with raw and corrected scores.

    Written by `export_pxd004242_pairs.py`, which checks that the run it
    exports reproduces the frozen rotation AUCs before saving.
    """
    path = root / "data/PXD004242/PXD004242_primary_pairs_rotation.tsv"
    table = pd.read_csv(path, sep="\t")
    labels = table["same_donor"].to_numpy(dtype=bool)

    def block(column: str) -> dict[str, object]:
        values = table[column].to_numpy(dtype=float)
        order = np.argsort(values)
        ranked = np.empty(len(order), dtype=float)
        ranked[order] = np.arange(1, len(order) + 1)
        pos, neg = labels, ~labels
        area = (ranked[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * neg.sum())
        return {"labels": labels, "values": values, "auc": float(area)}

    return {"raw": block("score_raw"), "corrected": block("score_corrected"),
            "n_features": 148}


def pxd004242(root: Path) -> dict[str, object]:
    base = root / "data/PXD004242"
    random273 = pd.read_csv(base / "PXD004242_random273_panel_rotation_summary.tsv", sep="\t")
    return {
        "adaptive": load_json(base / "PXD004242_random200_restricted.json"),
        "top273": load_json(base / "PXD004242_random200_top273.json"),
        "time": load_json(base / "PXD004242_time_distance_random200.json"),
        "swap": load_json(base / "PXD004242_swap_reference_curve_summary.json"),
        "random273": random273,
    }
