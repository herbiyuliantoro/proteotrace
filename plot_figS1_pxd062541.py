#!/usr/bin/env python3
"""Figure S1: the third cohort (PXD062541), run through the tool at its defaults.

A, B  pair scores before and after per-protein scaling, standardised against
      the different-donor distribution as in Figure 2. No batch information is
      deposited for this cohort, so there is no batch-correction step to show.
C     the coverage sensitivity analysis: AUC and same-person pairs found at
      the 99% cut-off for the default rule and three looser ones. Bars on the
      default are donor-cluster bootstrap 95% CIs.

Reads only tracked files: data/PXD062541_pairs.tsv and
data/PXD062541_primary.json, both written by analysis_pxd062541.py.
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

BINS = np.linspace(-4.0, 6.0, 41)
CENTRES = (BINS[:-1] + BINS[1:]) / 2


def auc(labels: np.ndarray, values: np.ndarray) -> float:
    order = np.argsort(values)
    ranked = np.empty(len(order), dtype=float)
    ranked[order] = np.arange(1, len(order) + 1)
    pos, neg = labels, ~labels
    return float((ranked[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * neg.sum()))


def load(root: Path):
    pairs = pd.read_csv(root / "data/PXD062541_pairs.tsv", sep="\t")
    summary = json.loads((root / "data/PXD062541_primary.json").read_text(encoding="utf-8"))
    labels = pairs["same_donor"].to_numpy(dtype=bool)
    frozen = summary["coverage"]["0.99"]
    for column, key in (("score_raw", "raw"), ("score_default", "default")):
        observed = auc(labels, pairs[column].to_numpy(dtype=float))
        if abs(observed - frozen[key]["auc"]) > 5e-4:
            raise SystemExit(f"{column}: AUC {observed:.4f} does not match the frozen "
                             f"{frozen[key]['auc']:.4f}; rerun analysis_pxd062541.py")
    return pairs, labels, summary


def draw_hist(ax, labels, raw_values, letter, title, show_legend):
    null = raw_values[~labels]
    values = (raw_values - null.mean()) / null.std()
    for mask, colour, name in ((~labels, S.GREY, "different people"),
                               (labels, S.AQUA, "same person")):
        y, _ = np.histogram(values[mask], bins=BINS, density=True)
        ax.fill_between(CENTRES, y, color=colour, alpha=0.35, linewidth=0, step="mid")
        ax.step(CENTRES, y, color=colour, linewidth=1.1, where="mid", label=name)
    ax.text(0.97, 0.80, f"AUC {auc(labels, raw_values):.3f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=9, fontweight="bold", color=S.INK)
    ax.set_xlim(BINS[0], BINS[-1])
    ax.set_xticks([-2, 0, 2, 4])
    ax.set_ylim(0, 1.35)          # headroom: legend and AUC sit above the bars
    ax.tick_params(axis="y", labelleft=False, left=False)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("similarity (standardised)")
    if show_legend:
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0))
    S.head(ax, letter, title)


def draw_coverage(ax, summary):
    keys = ["0.99", "0.95", "0.90", "0.80"]
    x = np.arange(len(keys))
    rows = [summary["coverage"][k] for k in keys]
    series = [("auc", S.AQUA, "AUC", -0.10), ("sens_99", S.ORANGE,
              "pairs found (99% cut-off)", 0.10)]
    handles = []
    for key, colour, label, dx in series:
        y = np.asarray([r["default"][key] for r in rows])
        ax.plot(x + dx, y, "o", color=colour, markersize=4.5, zorder=3)
        ci = rows[0]["default"]["ci"][key]
        ax.errorbar(x[0] + dx, y[0], yerr=[[y[0] - ci[0]], [ci[1] - y[0]]], fmt="none",
                    ecolor=colour, elinewidth=0.8, capsize=2, zorder=2)
        for xi, yi in zip(x + dx, y):
            ax.text(xi + 0.09 * np.sign(dx), yi, f"{yi:.2f}", ha="left" if dx > 0 else "right",
                    va="center", fontsize=7.5, color=S.INK)
        handles.append(Line2D([], [], marker="o", linestyle="none", color=colour,
                              markersize=4.5, label=label))
    labels = [f"≥{float(k):.0%}\n{r['features']}" for k, r in zip(keys, rows)]
    labels[0] = labels[0].replace("\n", "*\n")          # * = tool default (caption)
    ax.set_xticks(x, labels)
    ax.set_xlim(-0.6, len(keys) - 0.4)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("minimum share of runs a protein must be seen in")
    ax.set_ylabel("AUC  or  same-person pairs found")
    ax.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.70))
    S.head(ax, "C", "Using more proteins helps\nin this shallow cohort")


def make_figure(root: Path, out_dir: Path):
    pairs, labels, summary = load(root)
    fig = plt.figure(figsize=(S.WIDTH_2COL, 74 * S.MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.35], wspace=0.30,
                          left=0.075, right=0.975, top=0.83, bottom=0.21)
    draw_hist(fig.add_subplot(gs[0, 0]), labels, pairs["score_raw"].to_numpy(float),
              "A", "No correction", True)
    draw_hist(fig.add_subplot(gs[0, 1]), labels, pairs["score_default"].to_numpy(float),
              "B", "+ scale each protein\n(no batch information)", False)
    draw_coverage(fig.add_subplot(gs[0, 2]), summary)
    return S.save(fig, out_dir, "figS1_pxd062541")


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
