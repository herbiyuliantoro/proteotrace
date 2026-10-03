#!/usr/bin/env python3
"""Figure 5: the external cohort reproduces the result, and what it rests on.

A  every one of 200 random one-run assignments, through each preprocessing step
B  the feature panel, evaluated on the same 200 assignments: high-coverage
   proteins carry the signal; panels of arbitrary, mostly sparse proteins lose it
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


def draw_assignment(ax, draws, reference: float) -> None:
    modes = [("raw", "no\ncorrection", S.GREY), ("scaled", "scale each\nprotein", S.BLUE),
             ("corrected", "+ remove\nbatch effect", S.AQUA)]
    series = [draws[f"adaptive_{m}_primary_auc"].to_numpy(dtype=float) for m, _l, _c in modes]
    for a, b, c in zip(*series):
        ax.plot([0, 1, 2], [a, b, c], color=S.PALE, linewidth=0.5, alpha=0.25, zorder=1)
    parts = ax.violinplot(series, positions=[0, 1, 2], widths=0.5, showmeans=False,
                          showmedians=False, showextrema=False, bw_method=0.25)
    for body, (_m, _l, colour) in zip(parts["bodies"], modes):
        body.set_facecolor(colour)
        body.set_edgecolor("none")
        body.set_alpha(0.35)
    for x, values, (_m, _l, colour) in zip(range(3), series, modes):
        q025, med, q975 = np.quantile(values, [0.025, 0.5, 0.975])
        ax.vlines(x, q025, q975, color=S.INK, linewidth=0.7, zorder=4)
        ax.scatter(x, med, s=24, color=colour, edgecolor=S.SURFACE, linewidth=0.8, zorder=5)
        ax.text(x + 0.14, med, f"{med:.3f}", ha="left", va="center", fontsize=9,
                fontweight="bold", color=S.INK, zorder=6,
                bbox={"facecolor": S.SURFACE, "edgecolor": "none", "pad": 1.0})

    ax.axhline(reference, color=S.MUTED, linewidth=0.8, linestyle=(0, (3, 2)),
               zorder=2, label=f"twin cohort ({reference:.3f})")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, -0.01))
    ax.set_xticks([0, 1, 2], [l for _m, l, _c in modes])
    ax.set_xlim(-0.5, 2.75)
    ax.set_ylim(0.78, 1.0)
    ax.set_ylabel("AUC (1 = perfect)")
    S.head(ax, "A", "Same result for 200 random choices of run")


def draw_feature_depth(ax, summary: dict, adaptive_n: float) -> None:
    labels = [f"proteins seen\nin ~all runs ({adaptive_n:.0f})", "273 best-\nmeasured", "273 picked\nat random"]
    panels = ["adaptive", "top273", "random273"]
    x = np.arange(3)
    series = [("auc", S.AQUA, "AUC", -0.12),
              ("sens_99", S.ORANGE, "same-person pairs found (99% cut-off)", 0.12)]
    for key, colour, _label, dx in series:
        stats = [summary[f"{p}_corrected_primary_{key}"] for p in panels]
        med = np.asarray([s["median"] for s in stats])
        low = med - np.asarray([s["q025"] for s in stats])
        high = np.asarray([s["q975"] for s in stats]) - med
        ax.errorbar(x + dx, med, yerr=np.vstack([low, high]), fmt="none",
                    ecolor=colour, elinewidth=0.8, capsize=2, zorder=2)
        ax.scatter(x + dx, med, s=24, color=colour, zorder=3)
        for xi, yi in zip(x + dx, med):
            ax.text(xi + np.sign(dx) * 0.07, yi, f"{yi:.3f}",
                    ha="left" if dx > 0 else "right", va="center", fontsize=9,
                    color=S.INK, fontweight="bold")

    ax.set_xticks(x, labels)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(0.30, 1.04)
    ax.set_yticks(np.arange(0.3, 1.01, 0.1))
    ax.set_ylabel("AUC  or  same-person pairs found")
    handles = [Line2D([], [], marker="o", linestyle="none", color=c, markersize=4, label=lab)
               for _k, c, lab, _d in series]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.01))
    S.head(ax, "B", "Well-measured proteins work\nbetter than random ones")


def make_figure(root: Path, out_dir: Path):
    external = D.pxd004242(root)
    reference = D.load_json(root / "data/PXD001064_primary.json")[
        "estimands"]["primary_unrelated_cross_batch"]["corrected"]["auc"]
    fig = plt.figure(figsize=(S.WIDTH_2COL, 76 * S.MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1.1], wspace=0.26,
                          left=0.072, right=0.985, top=0.855, bottom=0.175)
    draw_assignment(fig.add_subplot(gs[0, 0]), external["draws"], reference)
    draw_feature_depth(fig.add_subplot(gs[0, 1]), external["primary"]["summary"],
                       external["primary"]["adaptive_features"]["median"])
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
