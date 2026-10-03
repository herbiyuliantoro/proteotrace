"""Claim 7: what the shipped tool's validation actually measures.

Uses the SHIPPED scoring/report code (that is the object of claim 7) but an
independent swap planter and independent accounting.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
import pandas as pd
from proteotrace.core import normalise, pair_table, similarity, threshold_from_pairs
from proteotrace.report import run_report

R = str(Path(__file__).resolve().parent.parent) + "/"
COH = {"PXD001064": ("data/PXD001064_protein_matrix.tsv", "examples/PXD001064_samples.tsv"),
       "PXD004242": ("data/PXD004242/PXD004242_protein_matrix.tsv", "examples/PXD004242_samples.tsv")}


def load(name):
    m, s = COH[name]
    X = pd.read_csv(R + m, sep="\t", index_col=0)
    S = pd.read_csv(R + s, sep="\t", dtype=str)
    S = S[S.run.isin(X.index)].reset_index(drop=True)
    return X.loc[S.run], S


def score_once(X, S):
    """The cohort is scored ONCE; labels only enter pair_table/threshold/report."""
    ann = S.set_index("run")
    nz = normalise(X, ann["batch"])
    return similarity(nz)


def verdicts(scores, S, margin, spec=0.99, rival_adjust=1):
    ann = S.set_index("run")
    pairs = pair_table(scores, ann["donor"], ann["batch"])
    match, rival = threshold_from_pairs(pairs, spec, rivals_per_run=rival_adjust)
    rep = run_report(scores, ann["donor"], ann["batch"], match, margin, rival_threshold=rival)
    return rep.set_index("run").loc[S.run]


# --- 1. what the repo's planter actually does
def repo_planter(S, n_pairs, rng):
    table, swapped = S.copy(), []
    for _, block in table.groupby("batch"):
        if len(swapped) >= n_pairs * 2:
            break
        idx = list(block.index); rng.shuffle(idx)
        for a, b in zip(idx[::2], idx[1::2]):
            if len(swapped) >= n_pairs * 2:
                break
            if table.loc[a, "donor"] != table.loc[b, "donor"]:
                table.loc[a, "donor"], table.loc[b, "donor"] = table.loc[b, "donor"], table.loc[a, "donor"]
                swapped += [table.loc[a, "run"], table.loc[b, "run"]]
    return table, swapped


def uniform_planter(S, n_pairs, rng):
    """Reciprocal within-batch swaps, pairs drawn uniformly over ALL batches."""
    table = S.copy()
    cands = []
    for _, block in table.groupby("batch"):
        idx = block.index.to_numpy()
        for a_i in range(len(idx)):
            for b_i in range(a_i + 1, len(idx)):
                a, b = idx[a_i], idx[b_i]
                if table.loc[a, "donor"] != table.loc[b, "donor"]:
                    cands.append((a, b))
    order = rng.permutation(len(cands))
    used, swapped = set(), []
    for k in order:
        a, b = cands[k]
        if a in used or b in used:
            continue
        da, db = table.loc[a, "donor"], table.loc[b, "donor"]
        # don't let two swaps touch the same donor, so collateral is attributable
        if any(table.loc[u, "donor"] in (da, db) for u in used):
            continue
        table.loc[a, "donor"], table.loc[b, "donor"] = db, da
        used |= {a, b}
        swapped += [table.loc[a, "run"], table.loc[b, "run"]]
        if len(swapped) >= 2 * n_pairs:
            break
    return table, swapped


def account(rep, S_true, S_planted, swapped):
    planted = rep.index.isin(swapped)
    flagged = (rep["verdict"] == "swap_suspected").to_numpy()
    assess = (rep["verdict"] != "no_reference").to_numpy()
    true_donor = S_true.set_index("run").loc[rep.index, "donor"]
    involved_donors = set(S_true.set_index("run").loc[swapped, "donor"])
    involved = true_donor.isin(involved_donors).to_numpy() & ~planted   # other runs of swapped donors
    clean = ~planted & ~involved
    f = lambda m: float((flagged & m & assess).sum() / max((m & assess).sum(), 1))
    return {"sens": f(planted), "fa_repo_def": f(~planted), "fa_collateral": f(involved),
            "fa_uninvolved": f(clean), "n_planted": int((planted & assess).sum()),
            "n_collateral": int((involved & assess).sum())}


if __name__ == "__main__":
    MARGINS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40]
    for name in COH:
        X, S = load(name)
        scores = score_once(X, S)
        print(f"\n######## {name}: runs={len(S)} donors={S.donor.nunique()} batches={S.batch.nunique()}")
        # (a) repo planter: which batches do the swaps come from?
        for rep_i in range(5):
            t, sw = repo_planter(S, 12, np.random.default_rng(100 + rep_i))
            b = S.set_index("run").loc[sw, "batch"].value_counts().to_dict()
            print(f"  repo planter seed {100+rep_i}: swapped runs by batch {b}")
        # (b) the clean cohort, nothing planted
        for m in (0.05, 0.20):
            v = verdicts(scores, S, m)["verdict"].value_counts().to_dict()
            print(f"  CLEAN cohort margin {m}: {v}")
        # (c) repo planter, decomposed
        rows = []
        for rep_i in range(5):
            t, sw = repo_planter(S, 12, np.random.default_rng(100 + rep_i))
            rows.append(account(verdicts(scores, t, 0.20), S, t, sw))
        print("  repo planter m=0.20 decomposed:", pd.DataFrame(rows).mean().round(4).to_dict())
        # (d) uniform planter, many repeats, margin sweep
        res = []
        for rep_i in range(40):
            rng = np.random.default_rng(5000 + rep_i)
            t, sw = uniform_planter(S, 12, rng)
            for m in MARGINS:
                a = account(verdicts(scores, t, m), S, t, sw)
                a.update(margin=m, rep=rep_i)
                res.append(a)
        R_ = pd.DataFrame(res)
        R_.to_csv(f"sweep_{name}.tsv", sep="\t", index=False)
        print("  uniform planter, 40 repeats x 12 pairs:")
        print(R_.groupby("margin")[["sens", "fa_repo_def", "fa_collateral", "fa_uninvolved"]].mean().round(4).to_string())
