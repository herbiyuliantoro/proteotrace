#!/usr/bin/env python3
"""Primary estimand, step ablation and donor-cluster CIs for both cohorts.

Written in response to the pre-preprint audit (docs/AUDIT_REPORT.md). Three
things changed relative to the earlier analysis scripts:

1. One estimand for both cohorts. Same-donor pairs on different batches
   against unrelated different-donor pairs on different batches. In PXD001064
   that excludes the 137 twin pairs (they are analysed separately in the
   relatedness ladder); in PXD004242 it keeps pairs across plate groups,
   because the tool scores rivals against every other batch. The earlier
   within-plate-group restriction is kept as a secondary, design-aware number.
2. A 2x2 ablation. The earlier "raw" baseline omitted both the per-batch
   median subtraction and the per-protein z-scaling, so the raw-to-corrected
   gain could not be attributed to either step. Every result is now reported
   for raw, scaling only, batch correction only, and both.
3. Uncertainty. Assignment quantiles describe run-selection variability only.
   Donor-cluster bootstrap CIs are added (twin pair is the cluster in
   PXD001064, donor in PXD004242).

    python analysis_primary_estimand.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from analysis_cross_batch import acquisition, read_metadata
from analysis_pxd004242_layer1 import read_sample_table, select_runs

TIME_WEEKS = {"-8W": -8, "0W": 0, "4W": 4, "13W": 13, "26W": 26, "39W": 39, "52W": 52}
MODES = {                      # (subtract per-batch median, z-scale each protein)
    "raw": (False, False),
    "scaled": (False, True),
    "batch": (True, False),
    "corrected": (True, True),
}
SPECIFICITIES = (0.95, 0.99)
DRAW_SEED = 20260919
RANDOM_PANEL_SEEDS = range(1000, 1020)


# ------------------------------------------------------------------ scoring
def log_centred(matrix: pd.DataFrame, features=None) -> pd.DataFrame:
    positive = matrix.replace(0, np.nan)
    if features is None:
        features = positive.columns[positive.notna().mean(axis=0) >= 0.99]
    values = np.log2(positive.loc[:, list(features)])
    return values.subtract(values.median(axis=1), axis=0)


def score(log_values: pd.DataFrame, batches: pd.Series, mode: str) -> np.ndarray:
    correct, scale = MODES[mode]
    x = log_values
    if correct:
        x = x - x.groupby(batches).transform("median")
    if scale:
        x = (x - x.mean(axis=0)) / x.std(axis=0).replace(0, np.nan)
    return x.T.corr(min_periods=30).to_numpy(dtype=float)


def operating(labels: np.ndarray, values: np.ndarray) -> dict[str, float]:
    negatives, positives = values[~labels], values[labels]
    out = {"auc": float(roc_auc_score(labels.astype(int), values))}
    for sp in SPECIFICITIES:
        out[f"sens_{int(sp * 100)}"] = float(np.mean(positives >= np.quantile(negatives, sp)))
    return out


def cluster_bootstrap(labels, values, left_unit, right_unit, n_boot=1000, seed=0):
    """Resample clusters with replacement; a pair is weighted by the product of
    its two clusters' multiplicities."""
    units = np.unique(np.r_[left_unit, right_unit])
    index = {u: k for k, u in enumerate(units)}
    li = np.asarray([index[u] for u in left_unit])
    ri = np.asarray([index[u] for u in right_unit])
    rng = np.random.default_rng(seed)
    out = {"auc": [], "sens_95": [], "sens_99": []}
    for _ in range(n_boot):
        w = rng.multinomial(len(units), np.full(len(units), 1 / len(units)))
        weight = w[li] * w[ri]
        keep = weight > 0
        take = np.repeat(np.flatnonzero(keep), weight[keep])
        y, v = labels[take], values[take]
        if y.all() or not y.any():
            continue
        for key, value in operating(y, v).items():
            out[key].append(value)
    return {k: [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))] for k, v in out.items()}


def quantiles(values) -> dict[str, float]:
    v = np.asarray([x for x in values if x is not None and np.isfinite(x)], dtype=float)
    return {"median": float(np.median(v)), "q025": float(np.quantile(v, 0.025)),
            "q975": float(np.quantile(v, 0.975))}


# ------------------------------------------------------------------ PXD001064
def pxd001064(root: Path):
    matrix = pd.read_csv(root / "data/PXD001064_protein_matrix.tsv", sep="\t", index_col=0)
    raw = read_metadata(root / "data/PXD001064_sample_info.txt")
    common = [run for run in matrix.index if run in raw]
    meta = pd.DataFrame([raw[run] for run in common], index=common)
    # one run per declared donor and visit; first in matrix row order
    meta = meta.loc[~meta.duplicated(["donor", "visit"], keep="first")].copy()
    meta["batch"] = [acquisition(run)[0] for run in meta.index]
    meta["pair_id"] = meta["donor"].str.extract(r"(Twin pair-\d+)-[AB]$")[0]
    return matrix.loc[meta.index], meta


def pxd001064_pairs(log_values, meta, mode):
    scores = score(log_values, meta["batch"], mode)
    ii, jj = np.triu_indices(len(meta), 1)
    v = scores[ii, jj]
    ok = np.isfinite(v)
    ii, jj, v = ii[ok], jj[ok], v[ok]
    donor, batch, pair = (meta[c].to_numpy() for c in ("donor", "batch", "pair_id"))
    zygosity = meta["zygosity"].to_numpy()
    same = donor[ii] == donor[jj]
    twin = (pair[ii] == pair[jj]) & ~same
    relation = np.where(same, "same_donor",
                        np.where(twin, np.where(zygosity[ii] == "MZ", "MZ_twins", "DZ_twins"),
                                 "unrelated"))
    return pd.DataFrame({"i": ii, "j": jj, "score": v, "same_donor": same,
                         "cross_batch": batch[ii] != batch[jj], "relation": relation,
                         "left_unit": pair[ii], "right_unit": pair[jj]})


def run_pxd001064(root: Path, n_boot: int) -> dict:
    matrix, meta = pxd001064(root)
    log_values = log_centred(matrix)
    result = {"runs": int(len(meta)), "donors": int(meta["donor"].nunique()),
              "features": int(log_values.shape[1]), "batches": int(meta["batch"].nunique()),
              "estimands": {}}
    for mode in MODES:
        pairs = pxd001064_pairs(log_values, meta, mode)
        cross = pairs.loc[pairs["cross_batch"]]
        primary = cross.loc[cross["relation"].isin(["same_donor", "unrelated"])]
        result["estimands"].setdefault("primary_unrelated_cross_batch", {})[mode] = operating(
            primary["same_donor"].to_numpy(), primary["score"].to_numpy())
        result["estimands"].setdefault("all_cross_batch_incl_twins", {})[mode] = operating(
            cross["same_donor"].to_numpy(), cross["score"].to_numpy())
        if mode == "corrected":
            p = primary
            result["n_positive"] = int(p["same_donor"].sum())
            result["n_negative"] = int((~p["same_donor"]).sum())
            result["ci_primary_corrected"] = cluster_bootstrap(
                p["same_donor"].to_numpy(), p["score"].to_numpy(),
                p["left_unit"].to_numpy(), p["right_unit"].to_numpy(), n_boot)
    return result


# ------------------------------------------------------------------ PXD004242
def pxd004242_pairs(scores: np.ndarray, meta: pd.DataFrame) -> pd.DataFrame:
    ii, jj = np.triu_indices(len(meta), 1)
    v = scores[ii, jj]
    ok = np.isfinite(v)
    ii, jj, v = ii[ok], jj[ok], v[ok]
    donor = meta["individual"].to_numpy()
    plate = meta["plate"].to_numpy()
    group = np.asarray([p[0] for p in plate])
    weeks = meta["sampling_time"].map(TIME_WEEKS).to_numpy()
    return pd.DataFrame({
        "i": ii, "j": jj, "score": v,
        "same_donor": donor[ii] == donor[jj],
        "cross_plate": plate[ii] != plate[jj],
        "same_group": group[ii] == group[jj],
        "gap_weeks": np.abs(weeks[ii] - weeks[jj]),
        "pre_to_post": (np.minimum(weeks[ii], weeks[jj]) < 0) & (np.maximum(weeks[ii], weeks[jj]) >= 0),
        "left_unit": donor[ii], "right_unit": donor[jj],
    })


def one_draw(matrix, rows, seed, top273, random_panel):
    meta = select_runs(rows, "random", None, seed)
    sub = matrix.loc[meta.index]
    record = {"seed": seed}
    for panel, features in (("adaptive", None), ("top273", top273), ("random273", random_panel)):
        log_values = log_centred(sub, features)
        if panel == "adaptive":
            record["adaptive_features"] = int(log_values.shape[1])
        for mode in MODES:
            pairs = pxd004242_pairs(score(log_values, meta["plate"], mode), meta)
            primary = pairs.loc[pairs["cross_plate"]]
            for key, value in operating(primary["same_donor"].to_numpy(),
                                        primary["score"].to_numpy()).items():
                record[f"{panel}_{mode}_primary_{key}"] = value
            within = primary.loc[primary["same_group"]]
            for key, value in operating(within["same_donor"].to_numpy(),
                                        within["score"].to_numpy()).items():
                record[f"{panel}_{mode}_withingroup_{key}"] = value
            if panel == "adaptive" and mode == "corrected":
                negatives = primary.loc[~primary["same_donor"], "score"].to_numpy()
                positives = primary.loc[primary["same_donor"]]
                record["n_positive"] = int(len(positives))
                record["n_negative"] = int(len(negatives))
                groups = list(positives.groupby("gap_weeks")) + [
                    ("pre_to_post", positives.loc[positives["pre_to_post"]]),
                    ("post_only", positives.loc[~positives["pre_to_post"]])]
                for gap, block in groups:
                    y = np.r_[np.ones(len(block), bool), np.zeros(len(negatives), bool)]
                    record[f"gap_{gap}_auc"] = float(roc_auc_score(y, np.r_[block["score"], negatives]))
                    record[f"gap_{gap}_n"] = int(len(block))
    return record


def run_pxd004242(root: Path, n_draws: int, n_boot: int) -> tuple[dict, pd.DataFrame]:
    base = root / "data/PXD004242"
    matrix = pd.read_csv(base / "PXD004242_protein_matrix.tsv", sep="\t", index_col=0)
    rows = read_sample_table(base / "PXD004242_biological_samples.tsv", "longitudinal")
    runs = [r["run"] for r in rows if r["run"] in matrix.index]
    coverage = matrix.loc[runs].replace(0, np.nan).notna().mean(axis=0)
    top273 = coverage.sort_values(ascending=False, kind="mergesort").index[:273].tolist()
    panels = [sorted(np.random.default_rng(s).choice(matrix.columns.to_numpy(), 273,
                                                     replace=False).tolist())
              for s in RANDOM_PANEL_SEEDS]
    records = []
    for k in range(n_draws):
        records.append(one_draw(matrix, rows, DRAW_SEED + k, top273, panels[k % len(panels)]))
        if (k + 1) % 20 == 0:
            print(f"completed_draws={k + 1}/{n_draws}", flush=True)
    draws = pd.DataFrame(records)

    key = "adaptive_corrected_primary_auc"
    target = draws[key].median()
    median_draw = int((draws[key] - target).abs().idxmin())
    median_seed = int(draws.loc[median_draw, "seed"])

    # donor-cluster CI on the median draw
    meta = select_runs(rows, "random", None, median_seed)
    pairs = pxd004242_pairs(score(log_centred(matrix.loc[meta.index]), meta["plate"], "corrected"), meta)
    primary = pairs.loc[pairs["cross_plate"]]
    ci = cluster_bootstrap(primary["same_donor"].to_numpy(), primary["score"].to_numpy(),
                           primary["left_unit"].to_numpy(), primary["right_unit"].to_numpy(), n_boot)

    summary = {
        "n_draws": int(n_draws),
        "draw_seed_start": DRAW_SEED,
        "random_panel_seeds": [RANDOM_PANEL_SEEDS.start, RANDOM_PANEL_SEEDS.stop - 1],
        "random_panel_rule": "draw k uses random panel seed 1000 + (k mod 20)",
        "median_draw": {"index": median_draw, "seed": median_seed,
                        **{c: float(draws.loc[median_draw, c]) for c in draws.columns
                           if c.startswith("adaptive_") and c.endswith(("_auc", "_sens_95", "_sens_99"))
                           and "_primary_" in c}},
        "ci_median_draw_corrected": ci,
        "adaptive_features": quantiles(draws["adaptive_features"]),
        "summary": {c: quantiles(draws[c]) for c in draws.columns
                    if c.endswith(("_auc", "_sens_95", "_sens_99"))},
        "gap_n_median": {c[4:-2]: float(draws[c].median()) for c in draws.columns
                         if c.startswith("gap_") and c.endswith("_n")},
        "interpretation": ("Quantiles describe run-assignment (and, for random273, panel) "
                           "variability, not confidence intervals; the donor-cluster bootstrap "
                           "CI on the median draw is the uncertainty estimate."),
    }
    return summary, draws


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument("--n-draws", type=int, default=200)
    ap.add_argument("--n-boot", type=int, default=1000)
    args = ap.parse_args()

    internal = run_pxd001064(args.root, args.n_boot)
    path = args.root / "data/PXD001064_primary.json"
    path.write_text(json.dumps(internal, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(internal["estimands"]["primary_unrelated_cross_batch"], indent=1))
    print("CI", internal["ci_primary_corrected"])

    external, draws = run_pxd004242(args.root, args.n_draws, args.n_boot)
    base = args.root / "data/PXD004242"
    (base / "PXD004242_primary_random200.json").write_text(
        json.dumps(external, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    draws.to_csv(base / "PXD004242_primary_random200.tsv", sep="\t", index=False)
    s = external["summary"]
    for mode in MODES:
        a = s[f"adaptive_{mode}_primary_auc"]
        print(f"PXD004242 {mode:9s} primary AUC median {a['median']:.4f} [{a['q025']:.4f}, {a['q975']:.4f}]"
              f"  sens99 {s[f'adaptive_{mode}_primary_sens_99']['median']:.4f}")
    print("median draw", external["median_draw"])
    print("CI", external["ci_median_draw_corrected"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
