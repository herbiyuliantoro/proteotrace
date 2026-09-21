#!/usr/bin/env python3
"""Figure 2: batch correction is what separates the two pair classes.

Raw and corrected scores live on different scales (raw correlations sit near
0.95 because overall abundance dominates), so a shared raw axis would hide the
comparison. Every panel is therefore standardised against its own
different-donor distribution: the null is centred at zero with unit spread, and
how far the same-donor distribution sits above it is what the method changes.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np

import figure_data as D

BINS = np.linspace(-4.0, 6.0, 51)
CENTRES = (BINS[:-1] + BINS[1:]) / 2


def draw(ax, block: dict, letter: str, title: str | None, cohort: str | None,
         note: str | None = None) -> None:
    labels, raw_values = block["labels"], block["values"]
    null = raw_values[~labels]
    values = (raw_values - null.mean()) / null.std()

    for mask, colour, name in ((~labels, S.GREY, "different donor"),
                               (labels, S.AQUA, "same donor")):
        y, _ = np.histogram(values[mask], bins=BINS, density=True)
        ax.fill_between(CENTRES, y, color=colour, alpha=0.35, linewidth=0, step="mid")
        ax.step(CENTRES, y, color=colour, linewidth=1.1, where="mid", label=name)

    gap = float(np.median(values[labels]) - np.median(values[~labels]))
    ax.text(0.97, 0.94, f"AUC {block['auc']:.4f}", transform=ax.transAxes, ha="right",
            va="top", fontsize=9, fontweight="bold", color=S.INK)
    ax.text(0.97, 0.80, f"gap {gap:.2f} SD", transform=ax.transAxes, ha="right",
            va="top", fontsize=9, color=S.MUTED)
    if note:
        ax.text(0.03, 0.94, note, transform=ax.transAxes, ha="left", va="top",
                fontsize=8, color=S.MUTED)

    ax.set_xlim(BINS[0], BINS[-1])
    ax.set_ylim(0, 0.62)
    ax.tick_params(axis="y", labelleft=False, left=False)
    if cohort:
        ax.set_ylabel(cohort, linespacing=1.6)
    S.ygrid(ax)
    S.head(ax, letter, title)


def make_figure(root: Path, out_dir: Path):
    internal = D.pxd001064_raw_vs_corrected(root)
    external = D.pxd004242_pair_scores(root)

    fig = plt.figure(figsize=(S.WIDTH_2COL, 112 * S.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.10, hspace=0.42,
                          left=0.135, right=0.985, top=0.905, bottom=0.095)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1], sharey=ax_a)
    ax_c = fig.add_subplot(gs[1, 0], sharey=ax_a)
    ax_d = fig.add_subplot(gs[1, 1], sharey=ax_a)

    draw(ax_a, internal["raw"], "A", "Before correction",
         "PXD001064 (DIA)\nall cross-batch pairs")
    draw(ax_b, internal["corrected"], "B", "After per-batch median subtraction", None)
    draw(ax_c, external["raw"], "C", None,
         "PXD004242 (DDA)\none rotation assignment")
    draw(ax_d, external["corrected"], "D", None, None,
         note="median over 200 assignments 0.9733")

    for ax in (ax_c, ax_d):
        ax.set_xlabel("score, in different-donor SD units")
    for ax in (ax_a, ax_b):
        ax.tick_params(axis="x", labelbottom=False)
    ax_a.legend(loc="upper left", bbox_to_anchor=(0.0, 0.99))
    return S.save(fig, out_dir, "fig2_batch_correction")


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
