#!/usr/bin/env python3
"""Figure 6: what the screen delivers in use.

A  donor signal against the gap between visits, including the pairs that span
   the weight-loss intervention (primary estimand, 200 assignments)
B  the shipped tool on planted swaps, against the number of visits each donor
   has: swaps caught, uninvolved runs flagged, and the swapped donors' other
   runs flagged. This replaces a closed-set, within-plate-group identification
   curve that was not what the tool does (audit F11).
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


def draw_time(ax, summary: dict, reference: float) -> None:
    gaps = sorted(int(k[4:-4]) for k in summary
                  if k.startswith("gap_") and k.endswith("_auc") and k[4:-4].isdigit())
    stats = [summary[f"gap_{g}_auc"] for g in gaps]
    x = np.asarray(gaps, dtype=float)
    y = np.asarray([s["median"] for s in stats])
    low = y - np.asarray([s["q025"] for s in stats])
    high = np.asarray([s["q975"] for s in stats]) - y
    ax.errorbar(x, y, yerr=np.vstack([low, high]), fmt="o", color=S.AQUA, ecolor=S.AQUA,
                linewidth=1.1, markersize=3.0, elinewidth=0.7, capsize=2, zorder=3)

    for key, position in (("pre_to_post", -14.0), ("post_only", -5.0)):
        s = summary[f"gap_{key}_auc"]
        med = s["median"]
        ax.errorbar(position, med, yerr=[[med - s["q025"]], [s["q975"] - med]],
                    fmt="o", color=S.INK, ecolor=S.INK, markersize=3.0,
                    elinewidth=0.7, capsize=2, zorder=3)
        ax.text(position, s["q025"] - 0.003, f"{med:.3f}", ha="center", va="top",
                fontsize=9, color=S.INK, fontweight="bold")
    ax.axvline(-1.0, color=S.GRID, linewidth=0.8, zorder=1)

    ax.axhline(reference, color=S.MUTED, linewidth=0.8, linestyle=(0, (3, 2)), zorder=2)
    ax.text(63, 0.902, f"Dashed line: PXD001064 ({reference:.3f})", ha="right",
            va="bottom", fontsize=7.5, color=S.MUTED)
    shown = [4, 13, 26, 39, 52, 60]
    ax.set_xticks([-14, -5] + shown, ["Across\ndiet", "After\ndiet"] + [str(v) for v in shown])
    ax.set_xlim(-19, 64)
    ax.set_ylim(0.90, 1.0)
    ax.set_xlabel("Interval between visits (weeks)")
    ax.set_ylabel("AUC")
    S.head(ax, "A", "Donor identity persists through the\nintervention and up to 60 weeks")


def draw_swap(ax, curve: dict) -> None:
    rows = curve["rows"]
    visits = np.asarray([row["n_visits"] for row in rows], dtype=float)
    # values at the line ends; identity from the legend, so the two
    # near-zero series do not need long labels that collide
    series = [("sensitivity", S.AQUA, "Swapped runs (sensitivity)", None),
              ("flagged_collateral", S.ORANGE, "Other runs of swapped donors", 0.10),
              ("flagged_uninvolved", S.BLUE, "Unaffected runs (false positives)", -0.005)]
    handles = []
    for key, colour, label, label_y in series:
        med = np.asarray([row[key]["mean"] for row in rows])
        low = med - np.asarray([row[key]["q025"] for row in rows])
        high = np.asarray([row[key]["q975"] for row in rows]) - med
        ax.errorbar(visits, med, yerr=np.vstack([low, high]), fmt="-o", color=colour,
                    ecolor=colour, markersize=3.0, linewidth=1.1, elinewidth=0.7,
                    capsize=2, zorder=3)
        ax.text(visits[-1] + 0.25, med[-1] if label_y is None else label_y,
                f"{med[-1]:.3f}", ha="left", va="center", fontsize=8.5, color=S.INK)
        handles.append(Line2D([], [], marker="o", color=colour, markersize=3.0,
                              linewidth=1.1, label=label))
    ax.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.55))

    ax.set_xlim(1.6, 8.0)
    ax.set_ylim(-0.05, 1.04)
    ax.set_xticks(np.arange(2, 8), [str(v) for v in range(2, 8)])
    ax.set_xlabel("Visits per donor")
    ax.set_ylabel("Fraction flagged as swap")
    S.head(ax, "B", "Swap detection improves with\nmore visits per donor")


def make_figure(root: Path, out_dir: Path):
    external = D.pxd004242(root)
    reference = D.load_json(root / "data/PXD001064_primary.json")[
        "estimands"]["primary_unrelated_cross_batch"]["corrected"]["auc"]
    fig = plt.figure(figsize=(S.WIDTH_2COL, 78 * S.MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.26,
                          left=0.072, right=0.985, top=0.855, bottom=0.17)
    draw_time(fig.add_subplot(gs[0, 0]), external["primary"]["summary"], reference)
    draw_swap(fig.add_subplot(gs[0, 1]), external["visit_curve"])
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
