#!/usr/bin/env python3
"""Figure 1: unsupervised batch correction and its measured effect.

A  analysis workflow
B  acquisition-day batch structure, PXD001064
C  plate batch structure, PXD004242
D  raw versus corrected AUC in both cohorts

Run from the ProteoTrace workspace:

    sh /mnt/c/Metabolomics_project/scripts/py plot_figure1_batch_correction.py
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from analysis_cross_batch import acquisition, read_metadata

BOX_FILL = "#f4f4f2"


def load_pxd001064(root: Path) -> tuple[dict[str, int], float, float]:
    matrix = pd.read_csv(root / "data/PXD001064_protein_matrix.tsv", sep="\t", index_col=0)
    raw_meta = read_metadata(root / "data/PXD001064_sample_info.txt")
    common = [run for run in matrix.index if run in raw_meta]
    matrix = matrix.loc[common]
    metadata = pd.DataFrame([raw_meta[run] for run in common], index=common)

    keep: list[str] = []
    seen: set[tuple[str, int]] = set()
    for run in matrix.index:
        key = (str(metadata.loc[run, "donor"]), int(metadata.loc[run, "visit"]))
        if key in seen:
            continue
        seen.add(key)
        keep.append(run)
    matrix = matrix.loc[keep]
    metadata = metadata.loc[keep].copy()
    metadata["batch"] = [acquisition(run)[0] for run in matrix.index]

    positive = matrix.replace(0, np.nan)
    selected = positive.notna().mean(axis=0) >= 0.99
    log_matrix = np.log2(positive.loc[:, selected])
    log_matrix = log_matrix.subtract(log_matrix.median(axis=1), axis=0)
    raw_scores = log_matrix.T.corr(min_periods=30)
    corrected = log_matrix - log_matrix.groupby(metadata["batch"]).transform("median")
    corrected = (corrected - corrected.mean(axis=0)) / corrected.std(axis=0).replace(0, np.nan)
    corrected_scores = corrected.T.corr(min_periods=30)

    def pair_auc(scores: pd.DataFrame) -> float:
        labels: list[int] = []
        values: list[float] = []
        names = list(scores.index)
        for i, left in enumerate(names):
            for right in names[i + 1:]:
                if metadata.loc[left, "batch"] == metadata.loc[right, "batch"]:
                    continue
                value = scores.loc[left, right]
                if np.isfinite(value):
                    labels.append(int(metadata.loc[left, "donor"] == metadata.loc[right, "donor"]))
                    values.append(float(value))
        order = np.argsort(values)
        ranked = np.empty(len(order), dtype=float)
        ranked[order] = np.arange(1, len(order) + 1)
        labels_arr = np.asarray(labels, dtype=int)
        pos = labels_arr == 1
        neg = ~pos
        return float((ranked[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * neg.sum()))

    return dict(sorted(Counter(metadata["batch"].tolist()).items())), pair_auc(raw_scores), pair_auc(corrected_scores)


def load_pxd004242(root: Path) -> tuple[dict[str, int], dict[str, float]]:
    samples = pd.read_csv(root / "data/PXD004242/PXD004242_biological_samples.tsv", sep="\t")
    longitudinal = samples.loc[samples["cohort"] == "longitudinal"]
    plate_counts: Counter[str] = Counter()
    for text in longitudinal["plates"]:
        plate_counts.update(str(text).split(","))
    payload = json.loads((root / "data/PXD004242/PXD004242_random200_restricted.json").read_text(encoding="utf-8"))
    s = payload["draw_summary"]
    aucs = {
        "raw": float(s["raw_primary_auc"]["median"]),
        "raw_q025": float(s["raw_primary_auc"]["q025"]),
        "raw_q975": float(s["raw_primary_auc"]["q975"]),
        "corrected": float(s["corrected_primary_auc"]["median"]),
        "corrected_q025": float(s["corrected_primary_auc"]["q025"]),
        "corrected_q975": float(s["corrected_primary_auc"]["q975"]),
    }
    return dict(sorted(plate_counts.items())), aucs


def draw_workflow(ax) -> None:
    ax.set_axis_off()
    S.head(ax, "A", letter_dx=-6.0)
    steps = [
        (0.045, "Run-level\nprotein matrix"),
        (0.385, "Sample-centre, then subtract\nper-batch protein medians"),
        (0.725, "Pairwise similarity and\ndonor-consistency score"),
    ]
    width, height, y = 0.23, 0.54, 0.34
    for x, text in steps:
        ax.add_patch(FancyBboxPatch(
            (x, y), width, height, boxstyle="round,pad=0.010,rounding_size=0.03",
            transform=ax.transAxes, linewidth=0.7, edgecolor="#b9b8b3", facecolor=BOX_FILL))
        ax.text(x + width / 2, y + height / 2, text, transform=ax.transAxes,
                ha="center", va="center", fontsize=6.5, color=S.INK, linespacing=1.35)
    for x in (0.285, 0.625):
        ax.add_patch(FancyArrowPatch(
            (x, y + height / 2), (x + 0.085, y + height / 2), transform=ax.transAxes,
            arrowstyle="-|>", mutation_scale=7, linewidth=0.8, color=S.MUTED))
    ax.text(0.045, 0.16, "PXD001064 DIA / OpenSWATH   ·   PXD004242 DDA / MaxQuant",
            transform=ax.transAxes, fontsize=6, color=S.MUTED, ha="left")
    ax.text(0.955, 0.16, "donor labels enter only at evaluation", transform=ax.transAxes,
            fontsize=6, color=S.ORANGE, ha="right")


def draw_day_batches(ax, counts: dict[str, int]) -> None:
    labels = list(counts)
    values = [counts[k] for k in labels]
    x = np.arange(len(labels))
    ax.bar(x, values, color=S.BLUE, width=0.68, linewidth=0)
    S.head(ax, "B", "PXD001064 acquisition-day batches",
           "231 runs · 116 donors · 273 protein groups")
    ax.set_ylabel("donor × visit runs")
    ax.set_xticks(x, [datetime.strptime(k, "%y%m%d").strftime("%d %b") for k in labels],
                  rotation=45, ha="right")
    ax.set_ylim(0, max(values) * 1.12)
    S.ygrid(ax)


def draw_plate_batches(ax, counts: dict[str, int]) -> None:
    cols, rows = list("ABCDE"), ["1", "2", "3", "4"]
    grid = np.full((len(rows), len(cols)), np.nan)
    for i, g in enumerate(rows):
        for j, s in enumerate(cols):
            if f"{g}{s}" in counts:
                grid[i, j] = counts[f"{g}{s}"]
    cmap = LinearSegmentedColormap.from_list("ptblue", list(S.SEQ_BLUE))
    cmap.set_bad(S.SURFACE)
    vmax = float(np.nanmax(grid))
    ax.imshow(grid, cmap=cmap, norm=Normalize(vmin=0, vmax=vmax), aspect="auto")
    S.head(ax, "C", "PXD004242 plate batches",
           "300 sample-timepoints · 4 runs per sample · 16 plates")
    ax.set_xticks(np.arange(len(cols)), cols)
    ax.set_yticks(np.arange(len(rows)), [f"group {r}" for r in rows])
    ax.set_xlabel("plate round")
    ax.set_ylabel("plate group")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if np.isfinite(grid[i, j]):
                shade = "white" if grid[i, j] / vmax > 0.55 else S.INK
                ax.text(j, i, f"{int(grid[i, j])}", ha="center", va="center", fontsize=6.5, color=shade)
            else:
                ax.text(j, i, "–", ha="center", va="center", fontsize=6.5, color=S.PALE)
    ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color=S.SURFACE, linewidth=1.4)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_auc(ax, raw_1: float, corrected_1: float, external: dict[str, float]) -> None:
    y = np.array([1.0, 0.0])
    raw = np.array([raw_1, external["raw"]])
    corrected = np.array([corrected_1, external["corrected"]])
    ax.set_xlim(0.74, 1.005)
    ax.set_ylim(-0.62, 1.62)
    for yi, rv, cv in zip(y, raw, corrected):
        ax.plot([rv, cv], [yi, yi], color=S.PALE, linewidth=2.4, solid_capstyle="round", zorder=1)
        ax.scatter(rv, yi, s=26, color=S.GREY, zorder=3)
        ax.scatter(cv, yi, s=30, color=S.AQUA, zorder=3)
        ax.text(rv, yi + 0.22, f"{rv:.3f}", ha="center", va="bottom", fontsize=6.5, color=S.MUTED)
        ax.text(cv, yi + 0.22, f"{cv:.3f}", ha="center", va="bottom", fontsize=6.5,
                color=S.AQUA, fontweight="bold")
        ax.text((rv + cv) / 2, yi - 0.20, f"+{cv - rv:.3f}", ha="center", va="top",
                fontsize=6.5, color=S.INK)
    ax.plot([external["raw_q025"], external["raw_q975"]], [0, 0], color=S.GREY, linewidth=1.4, zorder=2)
    ax.plot([external["corrected_q025"], external["corrected_q975"]], [0, 0], color=S.AQUA, linewidth=1.4, zorder=2)
    ax.set_yticks(y, ["PXD001064\nDIA / OpenSWATH", "PXD004242\nDDA / MaxQuant"])
    ax.set_xlabel("AUC, same-donor versus different-donor pairs")
    ax.set_xticks(np.arange(0.76, 1.01, 0.04))
    S.ygrid(ax, axis="x")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    S.head(ax, "D", "Unsupervised batch correction sets the operating signal",
           "PXD004242 bars: 2.5–97.5% run-assignment range across 200 assignments")
    handles = [Line2D([], [], marker="o", linestyle="none", color=S.GREY, markersize=4, label="raw"),
               Line2D([], [], marker="o", linestyle="none", color=S.AQUA, markersize=4,
                      label="batch-corrected")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.60, 1.00), ncol=2)


def make_figure(root: Path, out_dir: Path) -> tuple[Path, Path]:
    day_counts, raw_1, corrected_1 = load_pxd001064(root)
    plate_counts, external = load_pxd004242(root)

    fig = plt.figure(figsize=(S.WIDTH_2COL, 124 * S.MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.50, 1.30, 0.95], hspace=0.70, wspace=0.30,
                          left=0.135, right=0.985, top=0.935, bottom=0.085)
    draw_workflow(fig.add_subplot(gs[0, :]))
    draw_day_batches(fig.add_subplot(gs[1, 0]), day_counts)
    draw_plate_batches(fig.add_subplot(gs[1, 1]), plate_counts)
    draw_auc(fig.add_subplot(gs[2, :]), raw_1, corrected_1, external)
    return S.save(fig, out_dir, "figure1_batch_correction")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--out-dir", type=Path, default=Path("figures"))
    args = parser.parse_args()
    for path in make_figure(args.root, args.out_dir):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
