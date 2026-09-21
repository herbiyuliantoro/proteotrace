#!/usr/bin/env python3
"""Figure 4: what the method delivers as a quality-control screen.

A  sensitivity at fixed negative-pair specificity, both cohorts
B  within-plate swap screening against the number of reference visits held
   for each donor, which is what the headline swap rate depends on
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from analysis_cross_batch import pair_table
from analysis_mz_dz_ladder import prepare


def pxd001064_operating_points(root: Path) -> dict[str, float]:
    scores, metadata = prepare(root / "data/PXD001064_protein_matrix.tsv",
                               root / "data/PXD001064_sample_info.txt")
    pairs = pair_table(scores, metadata)
    pairs = pairs.loc[~pairs["same_batch"]]
    positives = pairs.loc[pairs["same_donor"], "score"].to_numpy()
    negatives = pairs.loc[~pairs["same_donor"], "score"].to_numpy()
    return {f"sens_{int(sp * 100)}": float(np.mean(positives >= float(np.quantile(negatives, sp))))
            for sp in (0.95, 0.99)}


def draw_operating_points(ax, internal: dict[str, float], adaptive: dict) -> None:
    x = np.asarray([0.0, 1.0])
    summary = adaptive["draw_summary"]
    internal_y = np.asarray([internal["sens_95"], internal["sens_99"]])
    external_y = np.asarray([summary["corrected_sens_95"]["median"],
                             summary["corrected_sens_99"]["median"]])
    external_lo = np.asarray([summary["corrected_sens_95"]["q025"],
                              summary["corrected_sens_99"]["q025"]])
    external_hi = np.asarray([summary["corrected_sens_95"]["q975"],
                              summary["corrected_sens_99"]["q975"]])

    ax.plot(x, internal_y, "-o", color=S.BLUE, markersize=3.6, linewidth=1.1, zorder=3)
    ax.errorbar(x, external_y, yerr=np.vstack([external_y - external_lo, external_hi - external_y]),
                fmt="-o", color=S.AQUA, ecolor=S.AQUA, markersize=3.6, linewidth=1.1,
                elinewidth=0.7, capsize=2, zorder=3)
    for xi, yi in zip(x, internal_y):
        dx, ha = (-0.06, "right") if xi == 0 else (0.06, "left")
        ax.text(xi + dx, yi - 0.012, f"{yi:.3f}", ha=ha, va="top", fontsize=6.5,
                color=S.BLUE, fontweight="bold")
    for xi, yi in zip(x, external_y):
        dx, ha = (-0.06, "right") if xi == 0 else (0.06, "left")
        ax.text(xi + dx, yi + 0.014, f"{yi:.3f}", ha=ha, va="bottom", fontsize=6.5,
                color=S.AQUA, fontweight="bold")

    ax.set_xticks(x, ["95% specificity", "99% specificity"])
    ax.set_xlim(-0.40, 1.40)
    ax.set_ylim(0.58, 0.96)
    ax.set_ylabel("same-donor sensitivity")
    S.ygrid(ax)
    S.head(ax, "A", "Sensitivity at a fixed false-flag rate",
           "PXD004242 bars: 2.5–97.5% run-assignment range")
    handles = [Line2D([], [], marker="o", color=S.BLUE, markersize=3.6, linewidth=1.1,
                      label="PXD001064 (DIA)"),
               Line2D([], [], marker="o", color=S.AQUA, markersize=3.6, linewidth=1.1,
                      label="PXD004242 (DDA)")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.01))


def draw_swap_curve(ax, curve: dict) -> None:
    rows = curve["rows"]
    visits = np.asarray([row["n_total_visits"] for row in rows], dtype=float)
    series = [
        ("swap_disagreement", S.AQUA, "flags the swap"),
        ("correct_swap_recovery", S.ORANGE, "names the true donor"),
    ]
    for stem, colour, label in series:
        med = np.asarray([row[f"{stem}_median"] for row in rows], dtype=float)
        lo = med - np.asarray([row[f"{stem}_q025"] for row in rows], dtype=float)
        hi = np.asarray([row[f"{stem}_q975"] for row in rows], dtype=float) - med
        ax.errorbar(visits, med, yerr=np.vstack([lo, hi]), fmt="-o", color=colour,
                    ecolor=colour, markersize=3.2, linewidth=1.1, elinewidth=0.7,
                    capsize=2, zorder=3)
        ax.text(visits[-1] + 0.18, med[-1], f"{label}\n{med[-1]:.3f}",
                ha="left", va="center", fontsize=6.3, color=colour, linespacing=1.35)
        ax.text(visits[0] - 0.12, med[0], f"{med[0]:.3f}  ", ha="right", va="center",
                fontsize=6.3, color=colour, fontweight="bold")

    ax.set_xlim(1.35, 9.9)
    ax.set_ylim(0.855, 1.012)
    ax.set_xticks(np.arange(2, 8), [str(v) for v in range(2, 8)])
    ax.set_xlabel("reference visits held for the donor")
    ax.set_ylabel("swap screening rate")
    S.ygrid(ax)
    S.head(ax, "B", "Swap screening depends on how many visits a cohort has",
           "PXD004242 · within-plate swaps · bars: 2.5–97.5% resampling range")


def make_figure(root: Path, out_dir: Path) -> tuple[Path, Path]:
    base = root / "data/PXD004242"
    adaptive = json.loads((base / "PXD004242_random200_restricted.json").read_text(encoding="utf-8"))
    curve = json.loads((base / "PXD004242_swap_reference_curve_summary.json").read_text(encoding="utf-8"))
    internal = pxd001064_operating_points(root)

    fig = plt.figure(figsize=(S.WIDTH_2COL, 74 * S.MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.55], wspace=0.24,
                          left=0.075, right=0.995, top=0.80, bottom=0.145)
    draw_operating_points(fig.add_subplot(gs[0, 0]), internal, adaptive)
    draw_swap_curve(fig.add_subplot(gs[0, 1]), curve)
    return S.save(fig, out_dir, "figure4_qc_operation")


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
