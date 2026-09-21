#!/usr/bin/env python3
"""Figure 6: what the screen delivers in use.

A  donor signal against the gap between visits, including the pair that spans
   the weight-loss intervention
B  swap screening against the number of reference visits a donor has, which is
   the condition the headline swap rate depends on
"""
from __future__ import annotations

import argparse
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np

import figure_data as D

REFERENCE = 0.9692


def draw_time(ax, time_summary: dict) -> None:
    rows = [row for row in time_summary["summary"] if row["category"] == "gap_weeks"]
    rows.sort(key=lambda row: int(row["gap"]))
    x = np.asarray([int(row["gap"]) for row in rows], dtype=float)
    y = np.asarray([row["auc_median"] for row in rows], dtype=float)
    low = y - np.asarray([row["auc_q025"] for row in rows], dtype=float)
    high = np.asarray([row["auc_q975"] for row in rows], dtype=float) - y
    ax.errorbar(x, y, yerr=np.vstack([low, high]), fmt="o-", color=S.AQUA, ecolor=S.AQUA,
                linewidth=1.1, markersize=3.0, elinewidth=0.7, capsize=2, zorder=3)

    special = {row["category"]: row for row in time_summary["summary"]
               if row["category"] in {"pre_to_post", "post_only"}}
    positions = {"pre_to_post": -11.0, "post_only": -4.5}
    for key, position in positions.items():
        row = special[key]
        med = row["auc_median"]
        ax.errorbar(position, med,
                    yerr=[[med - row["auc_q025"]], [row["auc_q975"] - med]],
                    fmt="o", color=S.INK, ecolor=S.INK, markersize=3.0,
                    elinewidth=0.7, capsize=2, zorder=3)
        ax.text(position, med + 0.006, f"{med:.3f}", ha="center", va="bottom",
                fontsize=9, color=S.INK, fontweight="bold")
    ax.axvline(-1.0, color=S.GRID, linewidth=0.8, zorder=1)

    ax.axhline(REFERENCE, color=S.GREY, linewidth=0.7, linestyle=(0, (3, 2)), zorder=2)
    shown = [4, 13, 26, 39, 52, 60]
    ax.set_xticks([-11, -4.5] + shown, ["pre→\npost", "post\nonly"] + [str(v) for v in shown])
    ax.set_xlim(-15, 64)
    ax.set_ylim(0.90, 1.0)
    ax.set_xlabel("between-visit gap (weeks)")
    ax.set_ylabel("cross-plate AUC")
    S.ygrid(ax)
    S.head(ax, "A", "Donor identity persists through the\nintervention and 60 weeks")


def draw_swap(ax, curve: dict) -> None:
    rows = curve["rows"]
    visits = np.asarray([row["n_total_visits"] for row in rows], dtype=float)
    series = [("swap_disagreement", S.AQUA, "label disagreement", "center", 0.0),
              ("correct_swap_recovery", S.ORANGE, "correct donor recovery", "center", 0.0)]
    for stem, colour, label, va, shift in series:
        med = np.asarray([row[f"{stem}_median"] for row in rows], dtype=float)
        low = med - np.asarray([row[f"{stem}_q025"] for row in rows], dtype=float)
        high = np.asarray([row[f"{stem}_q975"] for row in rows], dtype=float) - med
        ax.errorbar(visits, med, yerr=np.vstack([low, high]), fmt="-o", color=colour,
                    ecolor=colour, markersize=3.0, linewidth=1.1, elinewidth=0.7,
                    capsize=2, zorder=3)
        ax.text(visits[-1] + 0.25, med[-1] + shift, f"{label}\n{med[-1]:.3f}",
                ha="left", va=va, fontsize=9, color=colour, linespacing=1.4)

    ax.set_xlim(1.6, 13.0)
    ax.set_ylim(0.855, 1.035)
    ax.set_xticks(np.arange(2, 8), [str(v) for v in range(2, 8)])
    ax.set_xlabel("reference visits held for the donor")
    ax.set_ylabel("swap screening rate")
    S.ygrid(ax)
    S.head(ax, "B", "Simulated swaps: detection improves\nwith repeated sampling")


def make_figure(root: Path, out_dir: Path):
    external = D.pxd004242(root)
    fig = plt.figure(figsize=(S.WIDTH_2COL, 78 * S.MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.26,
                          left=0.072, right=0.985, top=0.855, bottom=0.17)
    draw_time(fig.add_subplot(gs[0, 0]), external["time"])
    draw_swap(fig.add_subplot(gs[0, 1]), external["swap"])
    return S.save(fig, out_dir, "fig6_qc_use")


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
