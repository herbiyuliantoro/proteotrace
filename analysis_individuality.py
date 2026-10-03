#!/usr/bin/env python3
"""Idea 2 of docs/EXTENSION_PLAN.md: weight proteins by how individual they are.

Each protein gets a weight equal to its intraclass correlation (ICC): the share
of its variance that lies between donors rather than within a donor. Two runs
are then compared with an ICC-weighted Pearson correlation instead of the plain
one. ICCs are always estimated on donors other than the ones whose pairs are
scored (design A), or on another cohort altogether (design B).

    python analysis_individuality.py            # development cohorts only
    python analysis_individuality.py --test     # adds the single held-out test
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

import analysis_primary_estimand as PE
from analysis_pxd004242_layer1 import read_sample_table, select_runs
from proteotrace.core import normalise
from proteotrace.io import read_matrix

MEDIAN_DRAW_SEED = 20261090
N_SPLITS = 20
N_BOOT = 1000


# ------------------------------------------------------------------ scoring
def weighted_corr(z: np.ndarray, w: np.ndarray, min_shared: int = 30) -> np.ndarray:
    """Weighted Pearson correlation between rows of z over pairwise-complete columns."""
    m = np.isfinite(z).astype(float)
    x = np.where(m > 0, z, 0.0)
    mw, xw = m * w, x * w
    sw = mw @ m.T
    sx = xw @ m.T
    sy = sx.T
    sxx = (x * xw) @ m.T
    syy = sxx.T
    sxy = xw @ x.T
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sxy - sx * sy / sw
        r = cov / np.sqrt((sxx - sx * sx / sw) * (syy - sy * sy / sw))
    r[(m @ m.T) < min_shared] = np.nan
    return r


def icc(z: np.ndarray, donors: np.ndarray) -> np.ndarray:
    """One-way random-effects ICC per column, from donors with >= 2 observations."""
    out = np.zeros(z.shape[1])
    for f in range(z.shape[1]):
        col = z[:, f]
        ok = np.isfinite(col)
        s = pd.Series(col[ok]).groupby(donors[ok])
        n_i = s.size()
        keep = n_i[n_i >= 2].index
        if len(keep) < 3:
            continue
        g = pd.Series(col[ok])[np.isin(donors[ok], keep)].groupby(donors[ok][np.isin(donors[ok], keep)])
        n_i, means = g.size(), g.mean()
        n, k = n_i.sum(), len(n_i)
        grand = g.sum().sum() / n
        ssb = (n_i * (means - grand) ** 2).sum()
        ssw = sum(((v - v.mean()) ** 2).sum() for _, v in g)
        msb, msw = ssb / (k - 1), ssw / (n - k)
        n0 = (n - (n_i ** 2).sum() / n) / (k - 1)
        sb = max((msb - msw) / n0, 0.0)
        out[f] = sb / (sb + msw) if (sb + msw) > 0 else 0.0
    return out


# --------------------------------------------------------------- evaluation
def auc(y, s):
    from scipy.stats import rankdata
    r = rankdata(s)
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def sens99(y, s):
    return float(np.mean(s[y] >= np.quantile(s[~y], 0.99)))


def metrics(y, s):
    return {"auc": auc(y, s), "sens_99": sens99(y, s)}


def paired_boot(y, s_new, s_old, left, right, n_boot=N_BOOT, seed=0):
    units = np.unique(np.r_[left, right])
    idx = {u: k for k, u in enumerate(units)}
    li = np.array([idx[u] for u in left]); ri = np.array([idx[u] for u in right])
    rng = np.random.default_rng(seed)
    d = {"auc": [], "sens_99": []}
    for _ in range(n_boot):
        c = rng.multinomial(len(units), np.full(len(units), 1 / len(units)))
        wt = c[li] * c[ri]
        take = np.repeat(np.flatnonzero(wt > 0), wt[wt > 0])
        yy = y[take]
        if yy.all() or not yy.any():
            continue
        a, b = metrics(yy, s_new[take]), metrics(yy, s_old[take])
        for k in d:
            d[k].append(a[k] - b[k])
    return {k: [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))] for k, v in d.items()}


# ------------------------------------------------------------------ cohorts
def first_accession(name: str) -> str | None:
    m = re.search(r"sp\|([A-Z0-9]+(?:-\d+)?)\|", name)
    if m:
        return m.group(1).split("-")[0]
    m = re.match(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})", name)
    return m.group(1) if m else None


def cohort_pxd001064(root):
    matrix, meta = PE.pxd001064(root)
    z = normalise(matrix, meta["batch"]).values
    donor = meta["donor"].to_numpy()
    unit = meta["pair_id"].to_numpy()
    batch = meta["batch"].to_numpy()
    ii, jj = np.triu_indices(len(meta), 1)
    same = donor[ii] == donor[jj]
    keep = (batch[ii] != batch[jj]) & (same | (unit[ii] != unit[jj]))     # unrelated negatives only
    return {"name": "PXD001064", "z": z, "donor": donor, "unit": unit, "batch": batch, "meta": meta, "ii": ii[keep], "jj": jj[keep],
            "same": same[keep], "acc": [first_accession(c) for c in z.columns]}


def cohort_pxd004242(root):
    base = root / "data/PXD004242"
    matrix = pd.read_csv(base / "PXD004242_protein_matrix.tsv", sep="\t", index_col=0)
    rows = read_sample_table(base / "PXD004242_biological_samples.tsv", "longitudinal")
    meta = select_runs(rows, "random", None, MEDIAN_DRAW_SEED)
    z = normalise(matrix.loc[meta.index], meta["plate"]).values
    donor = meta["individual"].astype(str).to_numpy()
    plate = meta["plate"].to_numpy()
    group = np.array([p[0] for p in plate])
    ii, jj = np.triu_indices(len(meta), 1)
    keep = plate[ii] != plate[jj]
    return {"name": "PXD004242", "z": z, "donor": donor, "unit": donor, "group": group, "batch": plate, "ii": ii[keep],
            "jj": jj[keep], "same": (donor[ii] == donor[jj])[keep],
            "acc": [first_accession(c.split(";")[0]) for c in z.columns]}


def cohort_pxd062541(root):
    matrix, _ = read_matrix(root / "data/PXD062541_report.tsv")
    samples = pd.read_csv(root / "examples/PXD062541_samples.tsv", sep="\t", dtype=str).set_index("run")
    matrix = matrix.loc[samples.index]
    z = normalise(matrix, samples["batch"]).values
    donor = samples["donor"].to_numpy()
    ii, jj = np.triu_indices(len(samples), 1)
    return {"name": "PXD062541", "z": z, "donor": donor, "unit": donor, "batch": samples["batch"].to_numpy(), "ii": ii, "jj": jj,
            "same": donor[ii] == donor[jj], "acc": [first_accession(c.split(";")[0]) for c in z.columns]}


# ------------------------------------------------------------------ designs
def evaluate(c, w):
    zz = c["z"].to_numpy(dtype=float)
    base = PE_corr(zz)[c["ii"], c["jj"]]
    new = weighted_corr(zz, w)[c["ii"], c["jj"]]
    ok = np.isfinite(base) & np.isfinite(new)
    return base[ok], new[ok], c["same"][ok], c["unit"][c["ii"][ok]], c["unit"][c["jj"][ok]]


def PE_corr(zz):
    return weighted_corr(zz, np.ones(zz.shape[1]))


def design_a(c, rng):
    zz = c["z"].to_numpy(dtype=float)
    units = np.unique(c["unit"])
    if "group" in c:                       # stratify by plate group
        halves = [set(), set()]
        for g in np.unique(c["group"]):
            u = np.unique(c["unit"][c["group"] == g]); rng.shuffle(u)
            halves[0] |= set(u[: len(u) // 2]); halves[1] |= set(u[len(u) // 2:])
    else:
        u = units.copy(); rng.shuffle(u)
        halves = [set(u[: len(u) // 2]), set(u[len(u) // 2:])]
    base_all, new_all, y_all, l_all, r_all = [], [], [], [], []
    base_full = PE_corr(zz)
    for train, test in ((halves[0], halves[1]), (halves[1], halves[0])):
        tr = np.isin(c["unit"], list(train))
        w = icc(zz[tr], c["donor"][tr])
        new_full = weighted_corr(zz, w)
        te = np.isin(c["unit"], list(test))
        sel = te[c["ii"]] & te[c["jj"]]
        b = base_full[c["ii"][sel], c["jj"][sel]]; n = new_full[c["ii"][sel], c["jj"][sel]]
        ok = np.isfinite(b) & np.isfinite(n)
        base_all.append(b[ok]); new_all.append(n[ok]); y_all.append(c["same"][sel][ok])
        l_all.append(c["unit"][c["ii"][sel][ok]]); r_all.append(c["unit"][c["jj"][sel][ok]])
    return tuple(np.concatenate(a) for a in (base_all, new_all, y_all, l_all, r_all))


def transfer_weights(source, target):
    zz = source["z"].to_numpy(dtype=float)
    w_src = icc(zz, source["donor"])
    by_acc = {}
    for a, wv in zip(source["acc"], w_src):
        if a and a not in by_acc:
            by_acc[a] = wv
    med = float(np.median(list(by_acc.values())))
    w = np.array([by_acc.get(a, med) for a in target["acc"]])
    matched = sum(a in by_acc for a in target["acc"])
    return w, matched


def report(label, base, new, y, l, r):
    mb, mn = metrics(y, base), metrics(y, new)
    ci = paired_boot(y, new, base, l, r)
    print(f"  {label}: baseline AUC {mb['auc']:.4f} sens99 {mb['sens_99']:.3f} | weighted AUC {mn['auc']:.4f} "
          f"sens99 {mn['sens_99']:.3f} | diff sens99 {mn['sens_99'] - mb['sens_99']:+.3f} CI {ci['sens_99'][0]:+.3f}..{ci['sens_99'][1]:+.3f}"
          f" | diff AUC {mn['auc'] - mb['auc']:+.4f} CI {ci['auc'][0]:+.4f}..{ci['auc'][1]:+.4f}  (pos {int(y.sum())}, neg {int((~y).sum())})")
    return {"baseline": mb, "weighted": mn, "diff_ci": ci, "n_pos": int(y.sum()), "n_neg": int((~y).sum())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument("--test", action="store_true", help="also run the single held-out test on PXD062541")
    ap.add_argument("--out", type=Path, default=Path("data/extension_individuality.json"))
    args = ap.parse_args()

    dev = [cohort_pxd001064(args.root), cohort_pxd004242(args.root)]
    result = {"plan": "docs/EXTENSION_PLAN.md", "development": {}, "test": None}
    for c in dev:
        print(f"\n{c['name']}: {c['z'].shape[0]} runs, {c['z'].shape[1]} proteins")
        b = PE_corr(c["z"].to_numpy(dtype=float))[c["ii"], c["jj"]]
        ok = np.isfinite(b)
        print(f"  sanity, unweighted primary AUC on all pairs: {auc(c['same'][ok], b[ok]):.4f}")
        diffs = []
        first = None
        for s in range(N_SPLITS):
            base, new, y, l, r = design_a(c, np.random.default_rng(1000 + s))
            diffs.append((metrics(y, new)["sens_99"] - metrics(y, base)["sens_99"],
                          metrics(y, new)["auc"] - metrics(y, base)["auc"]))
            if s == 0:
                first = report("design A, split 1", base, new, y, l, r)
        d = np.array(diffs)
        print(f"  design A over {N_SPLITS} splits: median diff sens99 {np.median(d[:, 0]):+.3f} "
              f"(range {d[:, 0].min():+.3f}..{d[:, 0].max():+.3f}); median diff AUC {np.median(d[:, 1]):+.4f}")
        result["development"][c["name"]] = {"design_A_split1": first,
                                            "design_A_median_diff_sens99": float(np.median(d[:, 0])),
                                            "design_A_median_diff_auc": float(np.median(d[:, 1])),
                                            "design_A_range_diff_sens99": [float(d[:, 0].min()), float(d[:, 0].max())]}
    w, matched = transfer_weights(dev[1], dev[0])
    print(f"\ndesign B: weights from PXD004242 -> PXD001064 ({matched}/{len(w)} proteins matched)")
    result["development"]["PXD001064"]["design_B"] = report("design B", *evaluate(dev[0], w))
    result["development"]["PXD001064"]["design_B"]["matched_proteins"] = matched

    if args.test:
        t = cohort_pxd062541(args.root)
        w, matched = transfer_weights(dev[1], t)
        print(f"\nHELD-OUT TEST design B: weights from PXD004242 -> PXD062541 ({matched}/{len(w)} matched)")
        result["test"] = report("PXD062541", *evaluate(t, w))
        result["test"]["matched_proteins"] = matched
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"\n{args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
