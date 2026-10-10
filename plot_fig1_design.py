#!/usr/bin/env python3
"""Figure 1: how the similarity score is computed, and which pairs it is read on.

A  the three steps, shown on a real slice of PXD001064: quantify, remove the
   per-batch offset, then correlate two runs
B  the pairing rule. Only pairs on different batches are read; across plate
   groups as well as within them, because the tool compares a run with every
   other batch. Same-batch pairs are never used
"""
from __future__ import annotations

import argparse
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import figure_data as D

DIVERGING = LinearSegmentedColormap.from_list("ptdiv", [S.BLUE, "#ffffff", S.ORANGE])


def draw_heatmap(ax, frame, batches, title: str, scale: float, show_ylabel: bool):
    values = frame.to_numpy(dtype=float).T          # proteins on rows
    ax.imshow(values, cmap=DIVERGING, norm=Normalize(-scale, scale), aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if show_ylabel:
        ax.set_ylabel("Proteins")

    # Mark the two acquisition batches under the columns.
    order = list(dict.fromkeys(batches))
    for index, batch in enumerate(order):
        columns = [i for i, b in enumerate(batches) if b == batch]
        left, right = min(columns) - 0.45, max(columns) + 0.45
        ax.plot([left, right], [values.shape[0] - 0.3] * 2, color=S.INK,
                linewidth=1.0, clip_on=False, transform=ax.transData)
        ax.text((left + right) / 2, values.shape[0] + 0.5, f"Batch {index + 1}",
                ha="center", va="top", fontsize=8, color=S.INK)
    ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left", pad=3)
    return values


def draw_scatter(ax, illustration):
    left, right = illustration["pair_scores"]
    mask = np.isfinite(left) & np.isfinite(right)
    ax.scatter(left[mask], right[mask], s=5, color=S.AQUA, alpha=0.55, linewidth=0)
    limit = 3.2
    ax.plot([-limit, limit], [-limit, limit], color=S.GREY, linewidth=0.7,
            linestyle=(0, (3, 3)), zorder=0)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xticks([-2, 0, 2])
    ax.set_yticks([-2, 0, 2])
    ax.set_xlabel("Run 1, batch 1 (z-score)")
    ax.set_ylabel("Run 2, batch 2 (z-score)")
    ax.set_title("3. Pairwise correlation", fontsize=8.5, fontweight="bold", loc="left", pad=3)
    ax.text(0.04, 0.96, f"Same donor\nr = {illustration['pair_r']:.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.5,
            color=S.INK, linespacing=1.4)


def draw_pairing(ax):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # One plate group shown in full, a second to show that pairs across groups count.
    ax.add_patch(FancyBboxPatch((0.045, 0.27), 0.46, 0.57,
                                boxstyle="round,pad=0.012,rounding_size=0.03",
                                linewidth=0.8, edgecolor="#b9b8b3", facecolor="#f7f7f5"))
    ax.add_patch(FancyBboxPatch((0.60, 0.27), 0.30, 0.57,
                                boxstyle="round,pad=0.012,rounding_size=0.03",
                                linewidth=0.8, edgecolor="#b9b8b3", facecolor="#f7f7f5"))
    ax.text(0.275, 0.895, "Plate group 1", ha="center", va="bottom", fontsize=8.5,
            color=S.INK, fontweight="bold")
    ax.text(0.750, 0.895, "Plate group 2", ha="center", va="bottom", fontsize=8.5,
            color=S.INK, fontweight="bold")

    columns = {"Plate B": 0.17, "Plate C": 0.385}
    rows = {"Donor 1": 0.715, "Donor 2": 0.455}
    for name, x in columns.items():
        ax.text(x, 0.375, name, ha="center", va="top", fontsize=8.5, color=S.MUTED)
    for name, y in rows.items():
        ax.text(0.055, y, name, ha="left", va="center", fontsize=8.5, color=S.INK)

    nodes = {}
    for row, y in rows.items():
        for column, x in columns.items():
            nodes[(row, column)] = (x, y)
            ax.scatter(x, y, s=58, facecolor=S.SURFACE, edgecolor=S.INK,
                       linewidth=0.9, zorder=3)
    other = (0.75, 0.585)
    ax.scatter(*other, s=58, facecolor=S.SURFACE, edgecolor=S.INK, linewidth=0.9, zorder=3)
    ax.text(0.75, 0.515, "Donor 3", ha="center", va="top", fontsize=8.5, color=S.INK)

    def link(a, b, colour, style, rad):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-", linewidth=1.3, color=colour,
                                     linestyle=style, shrinkA=5, shrinkB=5, zorder=2,
                                     connectionstyle=f"arc3,rad={rad}"))

    link(nodes[("Donor 1", "Plate B")], nodes[("Donor 1", "Plate C")], S.AQUA, "-", -0.30)
    link(nodes[("Donor 1", "Plate B")], nodes[("Donor 2", "Plate C")], S.GREY, "-", 0.22)
    link(nodes[("Donor 1", "Plate C")], other, S.GREY, "-", -0.12)
    link(nodes[("Donor 1", "Plate B")], nodes[("Donor 2", "Plate B")], S.MUTED, (0, (2, 2)), 0.35)

    handles = [
        Line2D([], [], color=S.AQUA, linewidth=1.3, label="Same donor, different plates (positive pair)"),
        Line2D([], [], color=S.GREY, linewidth=1.3,
               label="Different donors, different plates (negative pair)"),
        Line2D([], [], color=S.MUTED, linewidth=1.3, linestyle=(0, (2, 2)),
               label="Same plate (excluded)"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.10), ncol=1)
    S.head(ax, "B", "Pair definitions across plates", letter_dx=-22.0)


def make_figure(root: Path, out_dir: Path):
    illustration = D.pxd001064_illustration(root)
    raw, corrected = illustration["raw"], illustration["corrected"]
    scale = float(np.nanquantile(np.abs(raw.to_numpy(dtype=float)), 0.80))

    fig = plt.figure(figsize=(S.WIDTH_2COL, 118 * S.MM))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.05], hspace=0.78,
                             left=0.085, right=0.985, top=0.915, bottom=0.05)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.0, 1.0, 1.02], wspace=0.40)
    ax1 = fig.add_subplot(top[0, 0])
    ax2 = fig.add_subplot(top[0, 1])
    ax3 = fig.add_subplot(top[0, 2])

    draw_heatmap(ax1, raw, illustration["batches"], "1. Protein quantification", scale, True)
    draw_heatmap(ax2, corrected, illustration["batches"],
                 "2. Batch-median correction", scale, False)
    draw_scatter(ax3, illustration)
    S.head(ax1, "A", None, letter_dx=-26.0)

    # Short arrows just right of panels 1 and 2, in figure coordinates, so they
    # never run into the axis label of the next panel.
    for left_ax in (ax1, ax2):
        x0 = left_ax.get_position().x1
        y = (left_ax.get_position().y0 + left_ax.get_position().y1) / 2
        fig.patches.append(FancyArrowPatch(
            (x0 + 0.012, y), (x0 + 0.048, y), transform=fig.transFigure,
            arrowstyle="-|>", mutation_scale=12, linewidth=1.2, color=S.MUTED))

    box1, box2 = ax1.get_position(), ax2.get_position()
    cax = fig.add_axes((box1.x0 + 0.02, box1.y0 - 0.105, (box2.x1 - box1.x0) - 0.04, 0.013))
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-scale, scale), cmap=DIVERGING),
                       cax=cax, orientation="horizontal")
    bar.set_ticks([-scale, 0, scale])
    bar.set_ticklabels([f"-{scale:.1f}", "0", f"+{scale:.1f}"])
    bar.outline.set_visible(False)
    cax.tick_params(length=2, pad=1.5)
    cax.set_title("log2 intensity, median-centered (clipped)",
                  fontsize=8, fontweight="normal", color=S.INK, pad=3)

    draw_pairing(fig.add_subplot(outer[1]))
    return S.save(fig, out_dir, "fig1_design")


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
