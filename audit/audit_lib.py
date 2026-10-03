"""Independent re-implementation for the ProteoTrace audit.

Written from the method description, NOT by importing the repo's code.
Only raw inputs are read: the PXD001064 matrix + PRIDE sample map, the
PXD004242 MaxQuant zip + community SDRF.
"""
from __future__ import annotations

import io
import re
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from pathlib import Path

ROOT = str(Path(__file__).resolve().parent.parent)


# ------------------------------------------------------------------ metrics
def auc(y, s) -> float:
    """Mann-Whitney AUC with mid-ranks for ties."""
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    r = rankdata(s)
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def sens_at_spec(y, s, spec=0.99) -> float:
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    thr = np.quantile(s[~y], spec)
    return float(np.mean(s[y] >= thr))


# ------------------------------------------------------------------ scoring
def pairwise_corr(Z: np.ndarray, min_periods: int = 30) -> np.ndarray:
    """Pearson on pairwise-complete features, written with masked sums."""
    M = np.isfinite(Z).astype(float)
    X = np.where(M > 0, Z, 0.0)
    n = M @ M.T
    sx = X @ M.T                  # sum of x_i over features observed in both
    sy = sx.T
    sxx = (X * X) @ M.T
    syy = sxx.T
    sxy = X @ X.T
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sxy - sx * sy / n
        vx = sxx - sx * sx / n
        vy = syy - sy * sy / n
        r = cov / np.sqrt(vx * vy)
    r[n < min_periods] = np.nan
    return r


def preprocess(X: pd.DataFrame, batch: pd.Series, *, coverage=0.99, features=None,
               correct=True, fit=None, scale=True):
    """log2 -> run-median centre -> per-batch feature median -> z-score.

    `fit` is a boolean mask over rows: coverage, batch medians and z-scaling
    are estimated on those rows only and applied to all rows. Default: all.
    """
    if fit is None:
        fit = np.ones(len(X), bool)
    fit = np.asarray(fit, bool)
    P = X.where(X > 0)
    if features is None:
        cov = P.loc[fit].notna().mean(axis=0)
        features = cov.index[cov >= coverage]
    L = np.log2(P.loc[:, features].astype(float))
    L = L.sub(L.median(axis=1), axis=0)
    if correct:
        b = batch.loc[L.index]
        med = L.loc[fit].groupby(b.loc[fit]).median()
        # batches with no fit rows fall back to their own (unlabelled) median
        missing = sorted(set(b) - set(med.index))
        if missing:
            own = L.groupby(b).median().loc[missing]
            med = pd.concat([med, own])
        L = L - med.loc[b.to_numpy()].to_numpy()
    if scale:
        mu = L.loc[fit].mean(axis=0)
        sd = L.loc[fit].std(axis=0).replace(0, np.nan)
        L = (L - mu) / sd
    return L, list(features)


def score_matrix(X, batch, **kw):
    L, feats = preprocess(X, batch, **kw)
    return pairwise_corr(L.to_numpy()), feats


def pairs_from(S: np.ndarray, meta: pd.DataFrame, cols=("donor", "batch")):
    i, j = np.triu_indices(len(meta), 1)
    v = S[i, j]
    ok = np.isfinite(v)
    i, j, v = i[ok], j[ok], v[ok]
    out = {"i": i, "j": j, "score": v}
    for c in cols:
        a = meta[c].to_numpy()
        out["same_" + c] = a[i] == a[j]
    return pd.DataFrame(out)


# ------------------------------------------------------------------ PXD001064
def load_pxd001064(dedup="first"):
    text = open(f"{ROOT}/data/PXD001064_sample_info.txt", errors="replace").read()
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    rows = [l.split("\t") for l in lines[1:] if l.strip()]
    meta = pd.DataFrame(rows, columns=["run", "donor", "visit", "zygosity"])
    meta = meta[meta["visit"].isin(["1", "2"])].set_index("run")
    X = pd.read_csv(f"{ROOT}/data/PXD001064_protein_matrix.tsv", sep="\t", index_col=0)
    runs = [r for r in X.index if r in meta.index]
    meta = meta.loc[runs].copy()
    if dedup == "first":
        meta = meta[~meta.duplicated(["donor", "visit"], keep="first")]
    elif dedup == "last":
        meta = meta[~meta.duplicated(["donor", "visit"], keep="last")]
    meta["batch"] = [re.search(r"L(\d{6})_", r).group(1) for r in meta.index]
    meta["order"] = [int(re.search(r"L\d{6}_(\d+)", r).group(1)) for r in meta.index]
    meta["pair_id"] = meta["donor"].str.extract(r"(Twin pair-\d+)-[AB]")[0]
    return X.loc[meta.index], meta


# ------------------------------------------------------------------ PXD004242
def extract_pxd004242():
    with zipfile.ZipFile(f"{ROOT}/data/PXD004242/MaxQuant_Output.zip") as z:
        with z.open("proteinGroups.txt") as h:
            pg = pd.read_csv(h, sep="\t", low_memory=False, dtype=str)
    for flag in ("Only identified by site", "Reverse", "Potential contaminant"):
        pg = pg[pg[flag].fillna("") != "+"]
    ids = pg["Majority protein IDs"].where(pg["Majority protein IDs"].fillna("") != "", pg["Protein IDs"])
    cols = [c for c in pg.columns if c.startswith("Intensity ")]
    V = pg[cols].apply(pd.to_numeric)
    V.index = ids
    V.columns = [c[len("Intensity "):] for c in cols]
    return V.T


def load_pxd004242_meta():
    s = pd.read_csv(f"{ROOT}/data/PXD004242/PXD004242_community_annotated.sdrf.tsv",
                    sep="\t", dtype=str)
    f = s["comment[data file]"]
    m = f.str.extract(r"^(\d{8})_.*_Plate(\w+?)_(\d+)_(\d+)\.raw$")
    meta = pd.DataFrame({
        "run": m[1] + "_" + m[2] + "_" + m[3],
        "date": m[0],
        "plate": m[1],
        "donor": s["characteristics[individual]"],
        "time": s["characteristics[sampling time]"],
        "tech": s["comment[technical replicate]"],
        "label_donor": m[2],
        "label_tidx": m[3],
    })
    meta = meta[meta["run"].notna()].set_index("run")   # drops the 18 Top14 reference runs
    meta["group"] = meta["plate"].str[0]
    meta["sample"] = meta["donor"] + "|" + meta["time"]
    meta["batch"] = meta["plate"]
    return meta


TIME_ORDER = {"-8W": 1, "0W": 2, "4W": 3, "13W": 4, "26W": 5, "39W": 6, "52W": 7}
WEEKS = {"-8W": -8, "0W": 0, "4W": 4, "13W": 13, "26W": 26, "39W": 39, "52W": 52}


def longitudinal(meta):
    m = meta[meta["donor"].str.fullmatch(r"\d+")]
    n_t = m.groupby("donor")["time"].nunique()
    return m[m["donor"].isin(n_t.index[n_t > 2])].copy()


def draw_one_run_per_sample(meta, rng):
    pick = meta.groupby("sample", sort=True).apply(
        lambda b: b.index[rng.integers(0, len(b))], include_groups=False)
    return meta.loc[pick.to_numpy()]
