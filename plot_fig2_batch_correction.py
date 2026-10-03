#!/usr/bin/env python3
"""Figure 2: what each preprocessing step adds.

Three columns: raw (log2, run-median centred), plus per-protein scaling, plus
per-batch median subtraction. The earlier two-column version compared raw with
"after batch correction" while the second column also contained the scaling
step, which credited batch correction with gains that scaling makes (audit F1).

Raw and corrected scores live on different scales (raw correlations sit near
0.95 because overall abundance dominates), so a shared raw axis would hide the
comparison. Every panel is therefore standardised against its own
different-donor distribution: the null is centred at zero with unit spread, and
how far the same-donor distribution sits above it is what each step changes.
This is a per-panel affine map, so it leaves every AUC unchanged.
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


def draw(ax, block: dict, letter: str | None, title: str | None, cohort: str | None) -> None:
    labels, raw_values = block["labels"], block["values"]
    null = raw_values[~labels]
    values = (raw_values - null.mean()) / null.std()

    for mask, colour, name in ((~labels, S.GREY, "different people"),
                               (labels, S.AQUA, "same person")):
        y, _ = np.histogram(values[mask], bins=BINS, density=True)
        ax.fill_between(CENTRES, y, color=colour, alpha=0.35, linewidth=0, step="mid")
        ax.step(CENTRES, y, color=colour, linewidth=1.1, where="mid", label=name)

    ax.text(0.97, 0.94, f"AUC {block['auc']:.3f}", transform=ax.transAxes, ha="right",
            va="top", fontsize=9, fontweight="bold", color=S.INK)

    ax.set_xlim(BINS[0], BINS[-1])
    ax.set_xticks([-2, 0, 2, 4])            # inner ticks only: no clash between panels
    ax.set_ylim(0, 0.84)
    ax.tick_params(axis="y", labelleft=False, left=False)
    ax.spines["left"].set_visible(False)    # density has no scale worth reading
    if cohort:
        ax.set_ylabel(cohort, linespacing=1.6)
    if letter:
        S.head(ax, letter, title)


def make_figure(root: Path, out_dir: Path):
    internal = D.pxd001064_modes(root)
    external = D.pxd004242_pair_scores(root)

    fig = plt.figure(figsize=(S.WIDTH_2COL, 100 * S.MM))
    gs = fig.add_gridspec(2, 3, wspace=0.10, hspace=0.30,
                          left=0.095, right=0.99, top=0.89, bottom=0.10)
    axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]
    for row in axes:
        for ax in row[1:]:
            ax.sharey(row[0])

    columns = [("raw", "A", "No correction"), ("scaled", "B", "+ scale each protein"),
               ("corrected", "C", "+ remove batch effect")]
    cohorts = ["Twin cohort\n(PXD001064)", "Weight-loss cohort\n(PXD004242)"]
    for c, (mode, letter, title) in enumerate(columns):
        draw(axes[0][c], internal[mode], letter, title, cohorts[0] if c == 0 else None)
        draw(axes[1][c], external[mode], None, None, cohorts[1] if c == 0 else None)
        axes[1][c].set_xlabel("similarity (standardised)")
        axes[0][c].tick_params(axis="x", labelbottom=False)
    axes[0][0].legend(loc="upper left", bbox_to_anchor=(0.0, 0.99))
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
