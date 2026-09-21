#!/usr/bin/env python3
"""Figure 3: the external cohort reproduces the result, and what it depends on.

A  raw versus corrected AUC across 200 one-run assignments (PXD004242)
B  what the feature panel contributes: adaptive, top-coverage, random
C  donor signal against between-visit gap, including the intervention contrast
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import proteotrace_style as S

S.apply()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REFERENCE = 0.9692  # PXD001064 corrected cross-batch AUC


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def draw_assignment(ax, adaptive: dict) -> None:
    raw = np.asarray([d["raw"]["primary_within_group_cross_plate"]["auc"] for d in adaptive["draws"]], dtype=float)
    corrected = np.asarray([d["corrected"]["primary_within_group_cross_plate"]["auc"] for d in adaptive["draws"]], dtype=float)

    for a, b in zip(raw, corrected):
        ax.plot([0, 1], [a, b], color=S.PALE, linewidth=0.3, alpha=0.35, zorder=1)
    parts = ax.violinplot([raw, corrected], positions=[0, 1], widths=0.52, showmeans=False,
                          showmedians=False, showextrema=False, bw_method=0.25)
    for body, color in zip(parts["bodies"], [S.GREY, S.AQUA]):
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.35)
    for x, values, color in [(0, raw, S.GREY), (1, corrected, S.AQUA)]:
        q025, med, q975 = np.quantile(values, [0.025, 0.5, 0.975])
        ax.vlines(x, q025, q975, color=S.INK, linewidth=0.7, zorder=4)
        ax.scatter(x, med, s=26, color=color, edgecolor=S.SURFACE, linewidth=0.8, zorder=5)
        ax.text(x + 0.16, med, f"{med:.4f}", ha="left", va="center", fontsize=6.5,
                fontweight="bold", color=S.INK, zorder=6)

    ax.axhline(REFERENCE, color=S.ORANGE, linewidth=0.7, linestyle=(0, (3, 2)), zorder=2)
    ax.text(-0.52, REFERENCE + 0.004, f"PXD001064 reference {REFERENCE:.4f}", ha="left",
            va="bottom", fontsize=6, color=S.ORANGE)

    ax.set_xticks([0, 1], ["raw", "batch-corrected"])
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0.78, 1.0)
    ax.set_ylabel("cross-plate AUC")
    S.ygrid(ax)
    S.head(ax, "A", "External cohort: 200 run assignments",
           "one run per biological sample · lines pair one assignment")


def draw_feature_depth(ax, adaptive: dict, top273: dict, random_summary: pd.DataFrame) -> None:
    labels = ["adaptive\n(median 147)", "top-273\nby coverage", "random-273"]
    sources = [
        adaptive["draw_summary"],
        top273["draw_summary"],
        {
            "corrected_primary_auc": {
                "median": float(random_summary["corrected_auc"].median()),
                "q025": float(random_summary["corrected_auc"].quantile(0.025)),
                "q975": float(random_summary["corrected_auc"].quantile(0.975))},
            "corrected_sens_99": {
                "median": float(random_summary["sens99"].median()),
                "q025": float(random_summary["sens99"].quantile(0.025)),
                "q975": float(random_summary["sens99"].quantile(0.975))},
        },
    ]
    x = np.arange(3)
    series = [("corrected_primary_auc", S.AQUA, "AUC", -0.11),
              ("corrected_sens_99", S.ORANGE, "sensitivity at 99% specificity", 0.11)]
    for key, color, _label, dx in series:
        med = np.asarray([src[key]["median"] for src in sources], dtype=float)
        low = med - np.asarray([src[key]["q025"] for src in sources], dtype=float)
        high = np.asarray([src[key]["q975"] for src in sources], dtype=float) - med
        ax.errorbar(x + dx, med, yerr=np.vstack([low, high]), fmt="none",
                    ecolor=color, elinewidth=0.8, capsize=2, zorder=2)
        ax.scatter(x + dx, med, s=24, color=color, zorder=3)
        for xi, yi in zip(x + dx, med):
            ax.text(xi, yi + 0.028, f"{yi:.3f}", ha="center", va="bottom", fontsize=6.3,
                    color=color, fontweight="bold")

    ax.set_xticks(x, labels)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(0.40, 1.06)
    ax.set_yticks(np.arange(0.4, 1.01, 0.1))
    ax.set_ylabel("AUC or sensitivity")
    S.ygrid(ax)
    S.head(ax, "B", "Coverage quality matters more than panel size",
           "bars: 2.5–97.5% assignment or panel-selection range")
    handles = [Line2D([], [], marker="o", linestyle="none", color=c, markersize=4, label=lab)
               for _k, c, lab, _d in series]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.01), ncol=1)


def draw_time_distance(ax, time_summary: dict) -> None:
    rows = [row for row in time_summary["summary"] if row["category"] == "gap_weeks"]
    rows.sort(key=lambda row: int(row["gap"]))
    x = np.asarray([int(row["gap"]) for row in rows], dtype=float)
    y = np.asarray([row["auc_median"] for row in rows], dtype=float)
    low = y - np.asarray([row["auc_q025"] for row in rows], dtype=float)
    high = np.asarray([row["auc_q975"] for row in rows], dtype=float) - y
    ax.errorbar(x, y, yerr=np.vstack([low, high]), fmt="o-", color=S.AQUA, ecolor=S.AQUA,
                linewidth=1.1, markersize=3.2, elinewidth=0.7, capsize=2, zorder=3)

    special = {row["category"]: row for row in time_summary["summary"]
               if row["category"] in {"pre_to_post", "post_only"}}
    sx = {"pre_to_post": -9.0, "post_only": -4.0}
    for key, colour in (("pre_to_post", S.INK), ("post_only", S.INK)):
        row = special[key]
        med = row["auc_median"]
        ax.errorbar(sx[key], med,
                    yerr=[[med - row["auc_q025"]], [row["auc_q975"] - med]],
                    fmt="o", color=colour, ecolor=colour, markersize=3.2,
                    elinewidth=0.7, capsize=2, zorder=3)
        ax.text(sx[key], med + 0.006, f"{med:.3f}", ha="center", va="bottom",
                fontsize=6.3, color=S.INK, fontweight="bold")
    ax.axvline(-1.5, color=S.GRID, linewidth=0.8, zorder=1)

    ax.axhline(REFERENCE, color=S.GREY, linewidth=0.7, linestyle=(0, (3, 2)), zorder=2)
    ax.text(1.0, REFERENCE - 0.003, f"PXD001064 reference {REFERENCE:.4f}", ha="left",
            va="top", fontsize=6, color=S.MUTED)

    shown = [4, 13, 26, 39, 52, 60]
    ticks = [-9, -4] + shown
    names = ["pre→\npost", "post\nonly"] + [str(v) for v in shown]
    ax.set_xticks(ticks, names)
    ax.set_xlim(-12, 64)
    ax.set_ylim(0.90, 1.0)
    ax.set_xlabel("between-visit gap (weeks)")
    ax.set_ylabel("cross-plate AUC")
    ax.grid(axis="y", color=S.GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    S.head(ax, "C", "Donor signal survives the intervention and 60 weeks",
           "PXD004242 · error bars: 2.5–97.5% run-assignment range, not confidence intervals")


def make_figure(root: Path, out_dir: Path) -> tuple[Path, Path]:
    base = root / "data/PXD004242"
    adaptive = load_json(base / "PXD004242_random200_restricted.json")
    top273 = load_json(base / "PXD004242_random200_top273.json")
    time_summary = load_json(base / "PXD004242_time_distance_random200.json")
    random_summary = pd.read_csv(base / "PXD004242_random273_panel_rotation_summary.tsv", sep="\t")

    fig = plt.figure(figsize=(S.WIDTH_2COL, 118 * S.MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92], width_ratios=[1.0, 1.25],
                          hspace=0.62, wspace=0.26,
                          left=0.075, right=0.985, top=0.885, bottom=0.095)
    draw_assignment(fig.add_subplot(gs[0, 0]), adaptive)
    draw_feature_depth(fig.add_subplot(gs[0, 1]), adaptive, top273, random_summary)
    draw_time_distance(fig.add_subplot(gs[1, :]), time_summary)
    return S.save(fig, out_dir, "figure3_external_robustness")


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
