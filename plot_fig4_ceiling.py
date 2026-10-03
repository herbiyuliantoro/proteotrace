#!/usr/bin/env python3
"""Figure 4: where the discrimination stops, and why.

A, B  the two score distributions in each cohort (primary estimand), with the
      95% and 99% specificity thresholds drawn on them. The same-donor mass to
      the left of a threshold is what the screen misses.
C     the resulting operating points with donor-cluster bootstrap 95% CIs,
      which are much wider than the run-assignment ranges shown earlier.
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
    ax.set_xlabel("corrected similarity score")
    ax.tick_params(axis="y", labelleft=False, left=False)
    ax.set_ylabel("density")
    S.ygrid(ax)
    if show_legend:
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.87))
    S.head(ax, letter, title)
    return points


def draw_operating(ax, internal: dict, external: dict) -> None:
    x = np.asarray([0.0, 1.0])
    ci1 = internal["summary"]["ci_primary_corrected"]
    est1 = internal["summary"]["estimands"]["primary_unrelated_cross_batch"]["corrected"]
    summary = external["primary"]["summary"]
    ci2 = external["primary"]["ci_median_draw_corrected"]
    series = [
        ("PXD001064", S.BLUE, -0.05,
         np.asarray([est1["sens_95"], est1["sens_99"]]),
         np.asarray([ci1["sens_95"], ci1["sens_99"]])),
        ("PXD004242", S.AQUA, 0.05,
         np.asarray([summary["adaptive_corrected_primary_sens_95"]["median"],
                     summary["adaptive_corrected_primary_sens_99"]["median"]]),
         np.asarray([ci2["sens_95"], ci2["sens_99"]])),
    ]
    handles = []
    for label, colour, dx, y, ci in series:
        ax.errorbar(x + dx, y, yerr=np.vstack([y - ci[:, 0], ci[:, 1] - y]), fmt="-o",
                    color=colour, ecolor=colour, markersize=4.0, linewidth=1.2,
                    elinewidth=0.8, capsize=2, zorder=3)
        for xi, yi in zip(x + dx, y):
            outward = -1 if xi < 0.5 else 1          # 95% labels left, 99% labels right
            above = dx < 0                           # PXD001064 above, PXD004242 below
            ax.text(xi + outward * 0.10, yi + (0.012 if above else -0.012), f"{yi:.3f}",
                    ha="right" if outward < 0 else "left", va="bottom" if above else "top",
                    fontsize=9, color=S.INK, fontweight="bold")
        handles.append(Line2D([], [], marker="o", color=colour, markersize=4.0,
                              linewidth=1.2, label=label))

    ax.set_xticks(x, ["95%", "99%"])
    ax.set_xlim(-0.75, 1.65)
    ax.set_ylim(0.40, 1.0)
    ax.set_xlabel("negative-pair specificity")
    ax.set_ylabel("same-donor sensitivity")
    S.ygrid(ax)
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.01))
    S.head(ax, "C", "Operating points, 95% CI")


def make_figure(root: Path, out_dir: Path):
    internal = D.pxd001064_modes(root)
    external_pairs = D.pxd004242_pair_scores(root)
    external = D.pxd004242(root)

    fig = plt.figure(figsize=(S.WIDTH_2COL, 76 * S.MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 1.0], wspace=0.30,
                          left=0.075, right=0.985, top=0.855, bottom=0.165)
    draw_thresholds(fig.add_subplot(gs[0, 0]), internal["corrected"], "A",
                    "PXD001064 (DIA)", np.linspace(-0.55, 0.75, 53), True,
                    "unrelated cross-batch pairs")
    draw_thresholds(fig.add_subplot(gs[0, 1]), external_pairs["corrected"], "B",
                    "PXD004242 (DDA)", np.linspace(-0.55, 0.75, 53), False,
                    "median draw, cross-plate pairs")
    draw_operating(fig.add_subplot(gs[0, 2]), internal, external)
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
