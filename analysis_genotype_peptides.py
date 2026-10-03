#!/usr/bin/env python3
"""Idea 3 of docs/EXTENSION_PLAN.md: genotype-like peptides.

Uses the peptide-to-protein ratio matrices written by extract_peptide_ratios.py.
The peptide score is an ICC-weighted Pearson correlation over batch-corrected,
z-scaled ratios; the combined score is the mean of the current protein score and
the peptide score. ICC weights come from training donors (design A, cross-fitted)
or from PXD004242 (design B, matched by peptide sequence).

    python analysis_genotype_peptides.py               # development cohorts
    python analysis_genotype_peptides.py --test        # adds the held-out test
    python analysis_genotype_peptides.py --variants    # UniProt variant check
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from analysis_individuality import (N_SPLITS, cohort_pxd001064, cohort_pxd004242, cohort_pxd062541, icc, metrics,
                                    paired_boot, weighted_corr)

ROOT = Path(__file__).resolve().parent


def prep(ratio: pd.DataFrame, batch: np.ndarray) -> np.ndarray:
    x = ratio - ratio.groupby(batch).transform("median")
    x = (x - x.mean()) / x.std().replace(0, np.nan)
    return x.to_numpy(dtype=float)


def load_ratios(c):
    r = pd.read_csv(ROOT / f"data/peptides/{c['name']}_peptide_ratios.tsv.gz", sep="\t", index_col=0)
    runs = c["z"].index
    r = r.reindex(runs)
    return r, prep(r, np.asarray(c["batch"]))


def pair_scores(c, pep_matrix, w_pep, w_prot=None):
    zz = c["z"].to_numpy(dtype=float)
    prot = weighted_corr(zz, np.ones(zz.shape[1]))
    pep = weighted_corr(pep_matrix, w_pep, min_shared=30)
    out = {"protein": prot, "peptide": pep, "combined": (prot + pep) / 2}
    if w_prot is not None:
        out["combined_weighted_protein"] = (weighted_corr(zz, w_prot) + pep) / 2
    return out


def summarize(c, scores, sel=None, label=""):
    ii, jj, y = c["ii"], c["jj"], c["same"]
    if sel is not None:
        ii, jj, y = ii[sel], jj[sel], y[sel]
    base = scores["protein"][ii, jj]
    ok = np.isfinite(base)
    for k in scores:
        ok &= np.isfinite(scores[k][ii, jj])
    res = {"n_pos": int(y[ok].sum()), "n_neg": int((~y[ok]).sum()), "protein": metrics(y[ok], base[ok])}
    for k in scores:
        if k == "protein":
            continue
        s = scores[k][ii, jj][ok]
        res[k] = metrics(y[ok], s)
        res[k]["diff_ci"] = paired_boot(y[ok], s, base[ok], c["unit"][ii[ok]], c["unit"][jj[ok]])
        d = res[k]["sens_99"] - res["protein"]["sens_99"]
        ci = res[k]["diff_ci"]["sens_99"]
        print(f"  {label}{k:26s} AUC {res[k]['auc']:.4f} sens99 {res[k]['sens_99']:.3f} | vs protein "
              f"{res['protein']['sens_99']:.3f}: {d:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]")
    print(f"  {label}(pairs: {res['n_pos']} positive, {res['n_neg']} negative)")
    return res


def design_a(c, pep_matrix, rng):
    """Cross-fitted: ICC for peptides (and proteins) from one half, scored on the other."""
    units = np.unique(c["unit"])
    if "group" in c:
        halves = [set(), set()]
        for g in np.unique(c["group"]):
            u = np.unique(c["unit"][c["group"] == g]); rng.shuffle(u)
            halves[0] |= set(u[: len(u) // 2]); halves[1] |= set(u[len(u) // 2:])
    else:
        u = units.copy(); rng.shuffle(u)
        halves = [set(u[: len(u) // 2]), set(u[len(u) // 2:])]
    zz = c["z"].to_numpy(dtype=float)
    pooled = None
    sel_all = np.zeros(len(c["ii"]), bool)
    for train, test in ((halves[0], halves[1]), (halves[1], halves[0])):
        tr = np.isin(c["unit"], list(train))
        s = pair_scores(c, pep_matrix, icc(pep_matrix[tr], c["donor"][tr]), icc(zz[tr], c["donor"][tr]))
        te = np.isin(c["unit"], list(test))
        sel = te[c["ii"]] & te[c["jj"]]
        if pooled is None:
            pooled = {k: np.full_like(v, np.nan) for k, v in s.items()}
        for k, v in s.items():
            m = np.zeros_like(v, bool)
            m[c["ii"][sel], c["jj"][sel]] = True
            pooled[k][m] = v[m]
        sel_all |= sel
    return pooled, sel_all


def ladder(c, scores):
    meta = c["meta"]
    donor, pair, zyg, batch = (meta[k].to_numpy() for k in ("donor", "pair_id", "zygosity", "batch"))
    ii, jj = np.triu_indices(len(meta), 1)
    cross = batch[ii] != batch[jj]
    same = donor[ii] == donor[jj]
    twin = (pair[ii] == pair[jj]) & ~same
    rel = np.where(same, "same donor", np.where(twin, np.where(zyg[ii] == "MZ", "MZ twins", "DZ twins"), "unrelated"))
    out = {}
    for name in ("protein", "peptide", "combined"):
        s = scores[name][ii, jj]
        ok = cross & np.isfinite(s)
        thr = np.quantile(s[ok & (rel == "unrelated")], 0.99)
        out[name] = {r: {"median": float(np.median(s[ok & (rel == r)])),
                         "above_99": float(np.mean(s[ok & (rel == r)] >= thr))}
                     for r in ("same donor", "MZ twins", "DZ twins", "unrelated")}
        print(f"  ladder {name:9s} " + "  ".join(f"{r}: median {v['median']:+.3f}, above {v['above_99']:.1%}"
                                                 for r, v in out[name].items()))
    return out


def transfer(source_ratio_cols, w_source, target_cols):
    by_seq = dict(zip(source_ratio_cols, w_source))
    med = float(np.median(w_source))
    w = np.array([by_seq.get(s, med) for s in target_cols])
    return w, int(sum(s in by_seq for s in target_cols))


def variant_check(c42, pep42):
    """Do high-ICC peptides overlap known UniProt natural variants more often?"""
    info = pd.read_csv(ROOT / "data/peptides/PXD004242_peptide_info.tsv", sep="\t")
    w = icc(pep42, c42["donor"])
    info["icc"] = w
    cache = {}
    for acc in sorted(info["protein"].unique()):
        url = f"https://rest.uniprot.org/uniprotkb/{acc.split('-')[0]}.json?fields=sequence,ft_variant"
        try:
            d = json.load(urllib.request.urlopen(url, timeout=60))
        except Exception as e:  # noqa: BLE001
            print(f"  UniProt fetch failed for {acc}: {e}")
            continue
        seq = d.get("sequence", {}).get("value", "")
        var = []
        for f in d.get("features", []):
            if f.get("type") == "Natural variant":
                loc = f["location"]
                xref = [x["id"] for x in f.get("featureCrossReferences", []) if x.get("database") == "dbSNP"]
                var.append((loc["start"]["value"], loc["end"]["value"], f.get("description", ""), xref))
        cache[acc] = (seq, var)
        time.sleep(0.1)
    hits, snps = [], []
    for _, r in info.iterrows():
        seq, var = cache.get(r["protein"], ("", []))
        start = seq.find(r["sequence"]) + 1 if seq else 0
        if start <= 0:
            hits.append(np.nan); snps.append(""); continue
        end = start + len(r["sequence"]) - 1
        ov = [v for v in var if v[0] <= end and v[1] >= start]
        hits.append(bool(ov)); snps.append(";".join(x for v in ov for x in v[3]))
    info["overlaps_variant"] = hits
    info["dbsnp"] = snps
    known = info[info["overlaps_variant"].notna()].sort_values("icc", ascending=False)
    top = known.head(50)
    rest = known.iloc[50:]
    table = [[int(top["overlaps_variant"].sum()), int((~top["overlaps_variant"].astype(bool)).sum())],
             [int(rest["overlaps_variant"].sum()), int((~rest["overlaps_variant"].astype(bool)).sum())]]
    odds, p = fisher_exact(table, alternative="greater")
    print(f"  variant overlap: top 50 by ICC {table[0][0]}/50 ({table[0][0] / 50:.0%}); "
          f"other eligible peptides {table[1][0]}/{sum(table[1])} ({table[1][0] / sum(table[1]):.0%}); "
          f"Fisher one-sided p = {p:.2g}")
    known.head(15)[["sequence", "protein", "icc", "overlaps_variant", "dbsnp"]].to_csv(
        ROOT / "data/peptides/PXD004242_top_icc_peptides.tsv", sep="\t", index=False)
    return {"top50_overlap": table[0][0], "rest_overlap": table[1][0], "rest_n": sum(table[1]),
            "fisher_p_one_sided": float(p), "n_mapped": int(len(known))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--variants", action="store_true")
    ap.add_argument("--cohorts", default="PXD004242,PXD001064")
    ap.add_argument("--out", type=Path, default=ROOT / "data/extension_peptides.json")
    args = ap.parse_args()
    result = json.loads(args.out.read_text()) if args.out.exists() else {"development": {}, "test": None}

    c42 = cohort_pxd004242(ROOT)
    _, pep42 = load_ratios(c42)
    cols42 = pd.read_csv(ROOT / "data/peptides/PXD004242_peptide_ratios.tsv.gz", sep="\t", index_col=0, nrows=0).columns
    loaders = {"PXD004242": lambda: c42, "PXD001064": lambda: cohort_pxd001064(ROOT)}
    for name in args.cohorts.split(","):
        c = loaders[name]()
        ratio, pep = load_ratios(c)
        print(f"\n{name}: {pep.shape[1]} peptides; runs with peptide data {int(ratio.notna().any(axis=1).sum())}/{len(ratio)}")
        entry = {}
        diffs = {}
        for s in range(N_SPLITS):
            scores, sel = design_a(c, pep, np.random.default_rng(1000 + s))
            if s == 0:
                entry["design_A_split1"] = summarize(c, scores, sel, label="A split1 ")
            for k in ("peptide", "combined", "combined_weighted_protein"):
                ii, jj, y = c["ii"][sel], c["jj"][sel], c["same"][sel]
                b, v = scores["protein"][ii, jj], scores[k][ii, jj]
                ok = np.isfinite(b) & np.isfinite(v)
                diffs.setdefault(k, []).append(metrics(y[ok], v[ok])["sens_99"] - metrics(y[ok], b[ok])["sens_99"])
        for k, v in diffs.items():
            print(f"  design A, {N_SPLITS} splits, {k}: median diff sens99 {np.median(v):+.3f} "
                  f"(range {min(v):+.3f}..{max(v):+.3f})")
            entry[f"design_A_median_diff_sens99_{k}"] = float(np.median(v))
            entry[f"design_A_range_diff_sens99_{k}"] = [float(min(v)), float(max(v))]
        if name == "PXD001064":
            cols = ratio.columns
            w, matched = transfer(cols42, icc(pep42, c42["donor"]), cols)
            print(f"  design B: weights from PXD004242, {matched}/{len(cols)} peptides matched")
            scores = pair_scores(c, pep, w)
            entry["design_B"] = summarize(c, scores, label="B ")
            entry["design_B"]["matched_peptides"] = matched
            entry["ladder_design_B"] = ladder(c, scores)
        result["development"][name] = entry

    if args.variants:
        print("\nUniProt natural-variant check (PXD004242, all donors):")
        result["variant_check"] = variant_check(c42, pep42)

    if args.test:
        c = cohort_pxd062541(ROOT)
        ratio, pep = load_ratios(c)
        w, matched = transfer(cols42, icc(pep42, c42["donor"]), ratio.columns)
        print(f"\nHELD-OUT TEST PXD062541: {matched}/{len(ratio.columns)} peptides matched")
        result["test"] = summarize(c, pair_scores(c, pep, w), label="TEST ")
        result["test"]["matched_peptides"] = matched

    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"\n{args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
