#!/usr/bin/env python3
"""Figure 3: protein similarity by donor relationship.

Compare same-donor, MZ-twin, DZ-twin and unrelated pairs. The annotation marks
the overlap between same-donor and MZ-twin pairs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np

import figure_data as D

ORDER = ["same_donor", "MZ_twins", "DZ_twins", "unrelated"]
LABELS = ["Same\ndonor", "MZ twins\n(identical)", "DZ twins\n(fraternal)", "Unrelated"]
COLORS = [S.AQUA, S.BLUE, S.ORANGE, S.GREY]


def make_figure(root: Path, out_dir: Path):
    pairs = D.pxd001064_pairs(root)
    values = [pairs.loc[pairs["relation"] == rel, "score"].to_numpy() for rel in ORDER]
    x = np.arange(len(values))

    fig = plt.figure(figsize=(S.WIDTH_1COL, 78 * S.MM))      # JPR single column
    ax = fig.add_axes((0.175, 0.17, 0.80, 0.62))

    parts = ax.violinplot(values, positions=x, widths=0.76, showmeans=False,
                          showmedians=False, showextrema=False, bw_method=0.22)
    for body, colour in zip(parts["bodies"], COLORS):
        body.set_facecolor(colour)
        body.set_edgecolor("none")
        body.set_alpha(0.30)

    medians = np.array([float(np.median(v)) for v in values])
    for xi, arr, med, colour in zip(x, values, medians, COLORS):
        q05, q25, q75, q95 = np.quantile(arr, [0.05, 0.25, 0.75, 0.95])
        ax.vlines(xi, q05, q95, color=S.INK, linewidth=0.7, zorder=5)
        ax.vlines(xi, q25, q75, color=S.INK, linewidth=2.2, zorder=5)
        ax.scatter(xi, med, s=24, color=colour, edgecolor=S.SURFACE, linewidth=0.8, zorder=6)
        ax.text(xi + 0.13, med, f"{med:.2f}", ha="left", va="center", fontsize=7.5,
                color=S.INK, fontweight="bold", zorder=7)

    ax.axhline(0.0, color=S.GREY, linewidth=0.5, linestyle=(0, (1, 2)), zorder=0)   # r = 0 reference
    top = 0.70
    ax.plot([0, 0, 1, 1], [top - 0.03, top, top, top - 0.03], color=S.INK, linewidth=0.7)
    ax.text(0.5, top + 0.02, "MZ twins overlap with\nsame-donor pairs",
            ha="center", va="bottom", fontsize=7.5, color=S.INK, linespacing=1.35)

    ax.set_xticks(x, LABELS)
    ax.set_xlim(-0.6, 3.85)
    ax.set_ylim(-0.50, 0.90)
    ax.set_ylabel("Similarity (Pearson r)")
    S.head(ax, None, "Protein similarity by\ndonor relationship")
    return S.save(fig, out_dir, "fig3_relatedness")


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
