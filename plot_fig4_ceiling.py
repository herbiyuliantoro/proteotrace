#!/usr/bin/env python3
"""Figure 4: where the discrimination stops, and why.

A, B  the two score distributions in each cohort, with the 95% and 99%
      specificity thresholds drawn on them. The same-donor mass to the left of
      a threshold is what the screen misses, which is the honest reading of a
      sensitivity near 0.68.
C     the resulting operating points.
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


def draw_thresholds(ax, block: dict, letter: str, title: str, bins: np.ndarray,
                    show_legend: bool, estimand: str) -> dict:
    labels, values = block["labels"], block["values"]
    points = D.operating_points(labels, values)
    centres = (bins[:-1] + bins[1:]) / 2
    peak = 0.0
    for mask, colour, name in ((~labels, S.GREY, "different donor"),
                               (labels, S.AQUA, "same donor")):
        y, _ = np.histogram(values[mask], bins=bins, density=True)
        peak = max(peak, float(y.max()))
        ax.fill_between(centres, y, color=colour, alpha=0.32, linewidth=0, step="mid")
        ax.step(centres, y, color=colour, linewidth=1.1, where="mid", label=name)

    top = peak * 1.28
    for specificity, style, side in ((0.95, (0, (4, 2)), "right"), (0.99, "-", "left")):
        threshold = points[specificity]["threshold"]
        ax.axvline(threshold, color=S.ORANGE, linewidth=1.0, linestyle=style, zorder=4)
        pad = (bins[-1] - bins[0]) * (-0.012 if side == "right" else 0.012)
        ax.text(threshold + pad, top, f"{int(specificity * 100)}%", ha=side, va="top",
                fontsize=9, color=S.ORANGE, fontweight="bold")

    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, top * 1.02)
    ax.text(0.97, 0.04, estimand, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color=S.MUTED)
    ax.set_xlabel("corrected similarity score")
    ax.tick_params(axis="y", labelleft=False, left=False)
    ax.set_ylabel("density")
    S.ygrid(ax)
    if show_legend:
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.87))
    S.head(ax, letter, title)
    return points


def draw_operating(ax, internal: dict, adaptive: dict) -> None:
    x = np.asarray([0.0, 1.0])
    summary = adaptive["draw_summary"]
    internal_y = np.asarray([internal[0.95]["sensitivity"], internal[0.99]["sensitivity"]])
    external_y = np.asarray([summary["corrected_sens_95"]["median"],
                             summary["corrected_sens_99"]["median"]])
    lo = np.asarray([summary["corrected_sens_95"]["q025"], summary["corrected_sens_99"]["q025"]])
    hi = np.asarray([summary["corrected_sens_95"]["q975"], summary["corrected_sens_99"]["q975"]])

    ax.plot(x, internal_y, "-o", color=S.BLUE, markersize=4.0, linewidth=1.2, zorder=3)
    ax.errorbar(x, external_y, yerr=np.vstack([external_y - lo, hi - external_y]),
                fmt="-o", color=S.AQUA, ecolor=S.AQUA, markersize=4.0, linewidth=1.2,
                elinewidth=0.8, capsize=2, zorder=3)
    for values, colour, offset, va in ((internal_y, S.BLUE, -0.012, "top"),
                                       (external_y, S.AQUA, 0.014, "bottom")):
        for xi, yi in zip(x, values):
            dx, ha = (-0.07, "right") if xi == 0 else (0.07, "left")
            ax.text(xi + dx, yi + offset, f"{yi:.3f}", ha=ha, va=va, fontsize=9,
                    color=colour, fontweight="bold")

    ax.set_xticks(x, ["95%", "99%"])
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(0.58, 0.96)
    ax.set_xlabel("negative-pair specificity")
    ax.set_ylabel("same-donor sensitivity")
    S.ygrid(ax)
    handles = [Line2D([], [], marker="o", color=S.BLUE, markersize=4.0, linewidth=1.2,
                      label="PXD001064"),
               Line2D([], [], marker="o", color=S.AQUA, markersize=4.0, linewidth=1.2,
                      label="PXD004242, median of 200")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.01))
    S.head(ax, "C", "Operating points")


def make_figure(root: Path, out_dir: Path):
    internal = D.pxd001064_raw_vs_corrected(root)
    external_pairs = D.pxd004242_pair_scores(root)
    external = D.pxd004242(root)

    fig = plt.figure(figsize=(S.WIDTH_2COL, 76 * S.MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 1.0], wspace=0.30,
                          left=0.075, right=0.985, top=0.855, bottom=0.165)
    points = draw_thresholds(fig.add_subplot(gs[0, 0]), internal["corrected"], "A",
                             "PXD001064 (DIA)", np.linspace(-0.55, 0.75, 53), True,
                             "all cross-batch pairs")
    draw_thresholds(fig.add_subplot(gs[0, 1]), external_pairs["corrected"], "B",
                    "PXD004242 (DDA)", np.linspace(-0.55, 0.75, 53), False,
                    "one rotation assignment")
    draw_operating(fig.add_subplot(gs[0, 2]), points, external["adaptive"])
    return S.save(fig, out_dir, "fig4_ceiling")


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
