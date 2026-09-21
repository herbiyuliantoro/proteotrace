#!/usr/bin/env python3
"""Figure 5: the external cohort reproduces the result, and what it rests on.

A  every one of 200 arbitrary one-run assignments, raw against corrected
B  the feature panel: coverage quality, not panel size, carries the signal
"""
from __future__ import annotations

import argparse
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import figure_data as D

REFERENCE = 0.9692


def draw_assignment(ax, adaptive: dict) -> None:
    raw = np.asarray([d["raw"]["primary_within_group_cross_plate"]["auc"]
                      for d in adaptive["draws"]], dtype=float)
    corrected = np.asarray([d["corrected"]["primary_within_group_cross_plate"]["auc"]
                            for d in adaptive["draws"]], dtype=float)
    for a, b in zip(raw, corrected):
        ax.plot([0, 1], [a, b], color=S.PALE, linewidth=0.3, alpha=0.35, zorder=1)
    parts = ax.violinplot([raw, corrected], positions=[0, 1], widths=0.5, showmeans=False,
                          showmedians=False, showextrema=False, bw_method=0.25)
    for body, colour in zip(parts["bodies"], [S.GREY, S.AQUA]):
        body.set_facecolor(colour)
        body.set_edgecolor("none")
        body.set_alpha(0.35)
    for x, values, colour in ((0, raw, S.GREY), (1, corrected, S.AQUA)):
        q025, med, q975 = np.quantile(values, [0.025, 0.5, 0.975])
        ax.vlines(x, q025, q975, color=S.INK, linewidth=0.7, zorder=4)
        ax.scatter(x, med, s=24, color=colour, edgecolor=S.SURFACE, linewidth=0.8, zorder=5)
        ax.text(x + 0.14, med, f"{med:.4f}", ha="left", va="center", fontsize=9,
                fontweight="bold", color=S.INK, zorder=6)

    ax.axhline(REFERENCE, color=S.ORANGE, linewidth=0.9, linestyle=(0, (3, 2)),
               zorder=2, label="PXD001064")

    ax.legend(loc="lower right", bbox_to_anchor=(1.0, -0.01))
    ax.set_xticks([0, 1], ["raw", "batch-corrected"])
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0.78, 1.0)
    ax.set_ylabel("cross-plate AUC")
    S.ygrid(ax)
    S.head(ax, "A", "Stable across 200 run assignments")


def draw_feature_depth(ax, adaptive: dict, top273: dict, random273) -> None:
    labels = ["adaptive\n(median 147)", "top-273\nby coverage", "random-273"]
    sources = [
        adaptive["draw_summary"],
        top273["draw_summary"],
        {"corrected_primary_auc": {
            "median": float(random273["corrected_auc"].median()),
            "q025": float(random273["corrected_auc"].quantile(0.025)),
            "q975": float(random273["corrected_auc"].quantile(0.975))},
         "corrected_sens_99": {
            "median": float(random273["sens99"].median()),
            "q025": float(random273["sens99"].quantile(0.025)),
            "q975": float(random273["sens99"].quantile(0.975))}},
    ]
    x = np.arange(3)
    series = [("corrected_primary_auc", S.AQUA, "AUC", -0.10),
              ("corrected_sens_99", S.ORANGE, "sensitivity at 99% specificity", 0.10)]
    for key, colour, _label, dx in series:
        med = np.asarray([src[key]["median"] for src in sources], dtype=float)
        low = med - np.asarray([src[key]["q025"] for src in sources], dtype=float)
        high = np.asarray([src[key]["q975"] for src in sources], dtype=float) - med
        ax.errorbar(x + dx, med, yerr=np.vstack([low, high]), fmt="none",
                    ecolor=colour, elinewidth=0.8, capsize=2, zorder=2)
        ax.scatter(x + dx, med, s=24, color=colour, zorder=3)
        for xi, yi in zip(x + dx, med):
            ax.text(xi, yi + 0.028, f"{yi:.3f}", ha="center", va="bottom", fontsize=9,
                    color=colour, fontweight="bold")

    ax.set_xticks(x, labels)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(0.40, 1.06)
    ax.set_yticks(np.arange(0.4, 1.01, 0.1))
    ax.set_ylabel("AUC or sensitivity")
    S.ygrid(ax)
    handles = [Line2D([], [], marker="o", linestyle="none", color=c, markersize=4, label=lab)
               for _k, c, lab, _d in series]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.01))
    S.head(ax, "B", "High-coverage features outperform\nrandom panels of equal size")


def make_figure(root: Path, out_dir: Path):
    external = D.pxd004242(root)
    fig = plt.figure(figsize=(S.WIDTH_2COL, 76 * S.MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.26,
                          left=0.072, right=0.985, top=0.855, bottom=0.165)
    draw_assignment(fig.add_subplot(gs[0, 0]), external["adaptive"])
    draw_feature_depth(fig.add_subplot(gs[0, 1]), external["adaptive"],
                       external["top273"], external["random273"])
    return S.save(fig, out_dir, "fig5_transfer")


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
