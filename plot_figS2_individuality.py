#!/usr/bin/env python3
"""Figure S2: details of the individuality-weighted score.

A  change in sensitivity at 99% specificity when proteins are weighted by their
   intraclass correlation, with paired donor-cluster bootstrap 95% CIs
B  change in sensitivity at 99% specificity and in AUC across 20 random donor
   splits (cross-fitted weights) in the two development cohorts
C  per-protein intraclass correlation in PXD004242 against PXD001064 for the
   proteins quantified in both, which bounds how well weights transfer
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from matplotlib.lines import Line2D

from analysis_individuality import cohort_pxd001064, cohort_pxd004242, icc


def draw_weighting(ax, ind):
    dev = ind["development"]
    rows = [
        ("PXD004242\ncross-fitted weights", dev["PXD004242"]["design_A_split1"], S.AQUA),
        ("PXD001064\ncross-fitted weights", dev["PXD001064"]["design_A_split1"], S.AQUA),
        ("PXD001064\ntransferred weights", dev["PXD001064"]["design_B"], S.BLUE),
        ("PXD062541 (held-out)\ntransferred weights", ind["test"], S.ORANGE),
    ]
    y = np.arange(len(rows))[::-1]
    for yi, (label, r, colour) in zip(y, rows):
        d = r["weighted"]["sens_99"] - r["baseline"]["sens_99"]
        lo, hi = r["diff_ci"]["sens_99"]
        ax.errorbar(d, yi, xerr=[[d - lo], [hi - d]], fmt="o", color=colour, ecolor=colour,
                    markersize=5, elinewidth=1.0, capsize=2.5, zorder=3)
        ax.text(0.39, yi, f"{r['baseline']['sens_99']:.2f} → {r['weighted']['sens_99']:.2f}",
                ha="right", va="center", fontsize=8, color=S.INK)
    ax.text(0.39, len(rows) - 0.5, "Current → weighted", ha="right", va="center", fontsize=7.5, color=S.MUTED)
    ax.axvline(0, color=S.MUTED, linewidth=0.8, linestyle=(0, (3, 2)), zorder=1)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlim(-0.30, 0.40)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("Change in sensitivity at 99% specificity")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    S.head(ax, "A", "Change in sensitivity\nwith protein weighting", letter_dx=-106.5)



def draw_splits(ax, ind):
    dev = ind["development"]
    rows = [("PXD004242", "design_A_diffs_sens99", S.AQUA, 0), ("PXD001064", "design_A_diffs_sens99", S.AQUA, 1),
            ("PXD004242", "design_A_diffs_auc", S.GREY, 2), ("PXD001064", "design_A_diffs_auc", S.GREY, 3)]
    rng = np.random.default_rng(0)
    for name, key, colour, xi in rows:
        v = np.asarray(dev[name][key])
        ax.scatter(xi + rng.uniform(-0.12, 0.12, len(v)), v, s=12, color=colour, alpha=0.75,
                   linewidth=0, zorder=3)
        ax.hlines(np.median(v), xi - 0.22, xi + 0.22, color=S.INK, linewidth=1.2, zorder=4)
        ax.text(xi + 0.27, np.median(v), f"{np.median(v):+.3f}", va="center", fontsize=7.5, color=S.INK)
    ax.axhline(0, color=S.MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.set_xticks(range(4), ["PXD004242\nsensitivity", "PXD001064\nsensitivity",
                             "PXD004242\nAUC", "PXD001064\nAUC"])
    ax.set_xlim(-0.5, 3.8)
    ax.set_ylabel("Weighted minus current score")
    S.head(ax, "B", "Changes in AUC and sensitivity\nacross donor splits")


def draw_transfer(ax, root):
    c42, c64 = cohort_pxd004242(root), cohort_pxd001064(root)
    w42 = dict(zip(c42["acc"], icc(c42["z"].to_numpy(dtype=float), c42["donor"])))
    w64 = dict(zip(c64["acc"], icc(c64["z"].to_numpy(dtype=float), c64["donor"])))
    shared = [a for a in w42 if a and a in w64]
    x = np.array([w42[a] for a in shared]); y = np.array([w64[a] for a in shared])
    rho = spearmanr(x, y).statistic
    ax.scatter(x, y, s=12, color=S.BLUE, alpha=0.7, linewidth=0)
    ax.plot([0, 1], [0, 1], color=S.MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(0.97, 0.05, f"{len(shared)} shared proteins\nSpearman ρ = {rho:.2f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color=S.INK)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Protein ICC, PXD004242")
    ax.set_ylabel("Protein ICC, PXD001064")
    S.head(ax, "C", "Protein ICC across\ndevelopment cohorts")
    return len(shared), rho


def make_figure(root: Path, out_dir: Path):
    ind = json.loads((root / "data/extension_individuality.json").read_text())
    fig = plt.figure(figsize=(S.WIDTH_2COL, 140 * S.MM))
    # Panel A carries two-line row labels, so the top row gets a wider left margin;
    # its letter is offset to sit on the same vertical line as the letter of panel B.
    top = fig.add_gridspec(1, 1, left=0.24, right=0.975, top=0.9, bottom=0.6)
    bot = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.32, left=0.085, right=0.975, top=0.44, bottom=0.1)
    draw_weighting(fig.add_subplot(top[0, 0]), ind)
    draw_splits(fig.add_subplot(bot[0, 0]), ind)
    n, rho = draw_transfer(fig.add_subplot(bot[0, 1]), root)
    print(f"shared proteins {n}, spearman {rho:.3f}")
    return S.save(fig, out_dir, "figS2_individuality")


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
