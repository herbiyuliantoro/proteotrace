#!/usr/bin/env python3
"""Figure S3: genotype-associated peptide signatures.

A  overlap with UniProt natural-variant positions for the 50 peptides with the
   highest intraclass correlation in PXD004242, against all other peptides
B  median similarity by relatedness for the protein score and the peptide score,
   each relative to its same-donor median (PXD001064)
C  peptide-score distributions by relatedness in PXD001064 (weights transferred
   from PXD004242)
D  change in sensitivity at 99% specificity across 20 donor splits for the
   peptide score, the pre-registered combined score and the exploratory
   combination with the weighted protein score
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

from analysis_genotype_peptides import load_ratios, pair_scores, transfer
from analysis_individuality import cohort_pxd001064, cohort_pxd004242, icc

REL = ["same donor", "MZ twins", "DZ twins", "unrelated"]
REL_LABELS = ["Same\ndonor", "MZ\ntwins", "DZ\ntwins", "Unrelated"]
RELATIONS = REL
LABELS = ["Same\ndonor", "MZ twins\n(identical)", "DZ twins\n(fraternal)", "Unrelated"]
COLORS = [S.ORANGE] * 4  # one colour: every violin shows the peptide score, as in panel B


def draw_variants(ax, vc):
    top, rest_n, rest = vc["top50_overlap"], vc["rest_n"], vc["rest_overlap"]
    vals = [top / 50, rest / rest_n]
    ax.bar([0, 1], vals, width=0.55, color=[S.ORANGE, S.GREY], zorder=2)
    for xi, v, lab in zip([0, 1], vals, [f"{top}/50", f"{rest}/{rest_n}"]):
        ax.text(xi, v + 0.015, f"{v:.0%}\n({lab})", ha="center", va="bottom", fontsize=7.5, color=S.INK)
    mant, expo = f"{vc['fisher_p_one_sided']:.1e}".split("e")
    ax.text(0.97, 0.97, f"One-sided Fisher test\np = {mant} × 10$^{{{int(expo)}}}$", ha="right", va="top", fontsize=7.5,
            color=S.INK, transform=ax.transAxes)
    ax.set_xticks([0, 1], ["50 highest-ICC\npeptides", "Other eligible\npeptides"])
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, 0.75)
    ax.set_ylabel("Peptides overlapping a\nknown variant position")
    S.head(ax, "A", "Individual-specific peptides\nare enriched at variant sites")


def draw_ladder(ax, pep):
    lad = pep["development"]["PXD001064"]["ladder_design_B"]
    x = np.arange(len(RELATIONS))
    series = [("protein", "Protein score", S.BLUE, -0.12), ("peptide", "Peptide score", S.ORANGE, 0.12)]
    for key, label, colour, dx in series:
        ref = lad[key]["same donor"]["median"]
        vals = np.array([lad[key][r]["median"] / ref for r in RELATIONS])
        ax.plot(x + dx, vals, "o", color=colour, markersize=5.5, zorder=3)
        for xi, v in zip((x + dx)[1:], vals[1:]):  # the same-donor point is 1 by definition
            ax.text(xi + np.sign(dx) * 0.08, v, f"{v:.2f}".replace("-", "−"), ha="left" if dx > 0 else "right",
                    va="center", fontsize=7.5, color=S.INK)
    ax.axhline(1.0, color=S.MUTED, linewidth=0.8, linestyle=(0, (3, 2)), zorder=1)
    ax.axhline(0.0, color=S.GREY, linewidth=0.5, linestyle=(0, (1, 2)), zorder=1)
    ax.set_xticks(x, LABELS)
    ax.set_xlim(-0.6, len(RELATIONS) - 0.4)
    ax.set_ylim(-0.25, 1.6)
    ax.set_ylabel("Median similarity relative\nto same-donor pairs")
    handles = [Line2D([], [], marker="o", linestyle="none", color=c, markersize=5.5, label=l)
               for _k, l, c, _d in series]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 1.0))
    S.head(ax, "B", "Peptide similarity between\ntwins and repeat samples")



def draw_distributions(ax, root):
    c42, c64 = cohort_pxd004242(root), cohort_pxd001064(root)
    _, pep42 = load_ratios(c42)
    cols42 = pd.read_csv(root / "data/peptides/PXD004242_peptide_ratios.tsv.gz", sep="\t", index_col=0, nrows=0).columns
    ratio, pep = load_ratios(c64)
    w, _ = transfer(cols42, icc(pep42, c42["donor"]), ratio.columns)
    s = pair_scores(c64, pep, w)["peptide"]
    meta = c64["meta"]
    donor, pair, zyg, batch = (meta[k].to_numpy() for k in ("donor", "pair_id", "zygosity", "batch"))
    ii, jj = np.triu_indices(len(meta), 1)
    same = donor[ii] == donor[jj]; twin = (pair[ii] == pair[jj]) & ~same
    rel = np.where(same, "same donor", np.where(twin, np.where(zyg[ii] == "MZ", "MZ twins", "DZ twins"), "unrelated"))
    v = s[ii, jj]; ok = (batch[ii] != batch[jj]) & np.isfinite(v)
    data = [v[ok & (rel == r)] for r in REL]
    parts = ax.violinplot(data, positions=range(4), widths=0.75, showextrema=False, bw_method=0.25)
    for body, c in zip(parts["bodies"], COLORS):
        body.set_facecolor(c); body.set_edgecolor("none"); body.set_alpha(0.35)
    for xi, d in enumerate(data):
        q25, med, q75 = np.quantile(d, [0.25, 0.5, 0.75])
        ax.vlines(xi, q25, q75, color=S.INK, linewidth=2.0, zorder=4)
        ax.scatter(xi, med, s=18, color=COLORS[xi], edgecolor=S.SURFACE, linewidth=0.6, zorder=5)
        ax.text(xi + 0.12, med, f"{med:.2f}".replace("-", "−"), va="center", fontsize=7.5, color=S.INK)
    ax.axhline(0, color=S.GREY, linewidth=0.5, linestyle=(0, (1, 2)))
    ax.set_xticks(range(4), REL_LABELS)
    ax.set_xlim(-0.6, 3.75)
    ax.set_ylabel("Peptide score (Pearson r)")
    S.head(ax, "C", "Peptide scores by relatedness\n(PXD001064)")
    return {r: int(len(d)) for r, d in zip(REL, data)}


def draw_splits(ax, pep):
    keys = [("peptide", "Peptide\nscore"), ("combined", "Combined\n(prespecified)"),
            ("combined_weighted_protein", "Combined with\nweighted protein")]
    rng = np.random.default_rng(1)
    for ci, (name, colour) in enumerate((("PXD004242", S.AQUA), ("PXD001064", S.GREY))):
        for ki, (k, _lab) in enumerate(keys):
            v = np.asarray(pep["development"][name][f"design_A_diffs_sens99_{k}"])
            x = ki + (ci - 0.5) * 0.36
            ax.scatter(x + rng.uniform(-0.07, 0.07, len(v)), v, s=9, color=colour, alpha=0.75, linewidth=0, zorder=3,
                       label=name if ki == 0 else None)
            ax.hlines(np.median(v), x - 0.13, x + 0.13, color=S.INK, linewidth=1.1, zorder=4)
    ax.axhline(0, color=S.MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.set_xticks(range(3), [k[1] for k in keys])
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylabel("Change in sensitivity\nat 99% specificity")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.02))
    S.head(ax, "D", "Peptides do not improve\ndiscrimination")


def make_figure(root: Path, out_dir: Path):
    pep = json.loads((root / "data/extension_peptides.json").read_text())
    fig = plt.figure(figsize=(S.WIDTH_2COL, 150 * S.MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.85, 1.15], wspace=0.42, hspace=0.62,
                          left=0.1, right=0.985, top=0.9, bottom=0.1)
    draw_variants(fig.add_subplot(gs[0, 0]), pep["variant_check"])
    draw_ladder(fig.add_subplot(gs[0, 1]), pep)
    n = draw_distributions(fig.add_subplot(gs[1, 0]), root)
    draw_splits(fig.add_subplot(gs[1, 1]), pep)
    print("pairs by relation:", n)
    return S.save(fig, out_dir, "figS3_peptides")


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
