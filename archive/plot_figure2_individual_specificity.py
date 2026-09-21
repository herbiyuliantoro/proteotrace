#!/usr/bin/env python3
"""Figure 2: the similarity score tracks individual identity, not batch.

A  cross-batch score distributions with the relatedness ladder through the
   medians (same donor > MZ > DZ > unrelated)
B  same-donor ROC for PXD001064 with the 99%-specificity operating point
"""
from __future__ import annotations

import argparse
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, roc_curve

from analysis_cross_batch import pair_table
from analysis_mz_dz_ladder import prepare, relation

ORDER = ["same_donor", "MZ_twins", "DZ_twins", "unrelated"]
LABELS = ["same\ndonor", "MZ\ntwins", "DZ\ntwins", "unrelated"]
COLORS = [S.AQUA, S.BLUE, S.ORANGE, S.GREY]


def load_pairs(root: Path) -> pd.DataFrame:
    scores, metadata = prepare(
        root / "data/PXD001064_protein_matrix.tsv",
        root / "data/PXD001064_sample_info.txt",
    )
    pairs = pair_table(scores, metadata)
    pairs["relation"] = [relation(metadata.loc[left], metadata.loc[right])
                         for left, right in zip(pairs["left"], pairs["right"])]
    return pairs.loc[~pairs["same_batch"]].copy()


def draw_ladder(ax, pairs: pd.DataFrame) -> None:
    values = [pairs.loc[pairs["relation"] == rel, "score"].to_numpy() for rel in ORDER]
    x = np.arange(len(values))
    parts = ax.violinplot(values, positions=x, widths=0.74, showmeans=False,
                          showmedians=False, showextrema=False, bw_method=0.22)
    for body, color in zip(parts["bodies"], COLORS):
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.30)

    medians = np.array([float(np.median(v)) for v in values])
    ax.plot(x, medians, color=S.INK, linewidth=0.9, zorder=4)
    for xi, arr, med, color in zip(x, values, medians, COLORS):
        q05, q25, q75, q95 = np.quantile(arr, [0.05, 0.25, 0.75, 0.95])
        ax.vlines(xi, q05, q95, color=S.INK, linewidth=0.7, zorder=5)
        ax.vlines(xi, q25, q75, color=S.INK, linewidth=2.2, zorder=5)
        ax.scatter(xi, med, s=26, color=color, edgecolor=S.SURFACE, linewidth=0.8, zorder=6)
        ax.text(xi + 0.16, med, f"{med:.3f}", ha="left", va="center", fontsize=6.5,
                color=S.INK, fontweight="bold", zorder=7)

    ax.axhline(0.0, color=S.GRID, linewidth=0.6, zorder=0)
    ax.set_xticks(x, [f"{lab}\nn={len(v):,}" for lab, v in zip(LABELS, values)])
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(-0.50, 0.80)
    ax.set_ylabel("cross-batch similarity score")
    S.ygrid(ax)
    S.head(ax, "A", "Score tracks relatedness, not batch",
           "PXD001064 · cross-batch pairs only · bars: 5–95% and 25–75%")

    # The overlap between same-donor and MZ pairs is the tool's blind spot, so
    # it is marked on the panel rather than left to the caption.
    top = 0.70
    ax.plot([0, 0, 1, 1], [top - 0.03, top, top, top - 0.03], color=S.ORANGE, linewidth=0.7)
    ax.text(0.5, top + 0.015, "overlapping: identical twins are not separable",
            ha="center", va="bottom", fontsize=6, color=S.ORANGE)


def draw_roc(ax, pairs: pd.DataFrame) -> None:
    y_true = pairs["same_donor"].astype(int).to_numpy()
    score = pairs["score"].to_numpy()
    fpr, tpr, _ = roc_curve(y_true, score)
    roc_auc = auc(fpr, tpr)

    ax.plot([0, 1], [0, 1], color=S.PALE, linewidth=0.8, linestyle=(0, (3, 3)), zorder=1)
    ax.plot(fpr, tpr, color=S.AQUA, linewidth=1.6, zorder=3)

    different = pairs.loc[~pairs["same_donor"], "score"].to_numpy()
    same = pairs.loc[pairs["same_donor"], "score"].to_numpy()
    threshold = float(np.quantile(different, 0.99))
    sensitivity = float(np.mean(same >= threshold))
    fpr_op = float(np.mean(different >= threshold))
    ax.scatter(fpr_op, sensitivity, s=26, color=S.ORANGE, zorder=5)
    ax.annotate(f"99% specificity\nsensitivity {sensitivity:.3f}",
                xy=(fpr_op, sensitivity), xytext=(fpr_op + 0.07, sensitivity - 0.10),
                fontsize=6.5, color=S.ORANGE, ha="left", va="top",
                arrowprops={"arrowstyle": "-", "color": S.ORANGE, "lw": 0.6,
                            "shrinkA": 1, "shrinkB": 2})

    ax.text(0.97, 0.06, f"AUC = {roc_auc:.4f}", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=8, fontweight="bold", color=S.AQUA)
    ax.text(0.80, 0.76, "chance", ha="left", va="bottom", fontsize=6.5, color=S.MUTED)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.01)
    ax.set_xticks(np.arange(0, 1.01, 0.25))
    ax.set_yticks(np.arange(0, 1.01, 0.25))
    ax.set_xlabel("false-positive rate")
    ax.set_ylabel("true-positive rate")
    ax.grid(color=S.GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    S.head(ax, "B", "Same-donor discrimination, PXD001064",
           f"{int(y_true.sum()):,} same-donor and {int((y_true == 0).sum()):,} different-donor pairs")


def make_figure(root: Path, out_dir: Path) -> tuple[Path, Path]:
    pairs = load_pairs(root)
    fig = plt.figure(figsize=(S.WIDTH_2COL, 78 * S.MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.28,
                          left=0.075, right=0.985, top=0.80, bottom=0.145)
    draw_ladder(fig.add_subplot(gs[0, 0]), pairs)
    draw_roc(fig.add_subplot(gs[0, 1]), pairs)
    return S.save(fig, out_dir, "figure2_individual_specificity")


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
