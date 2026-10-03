"""Claims 5 (time gap) and 6 (panel quality), re-derived."""
import numpy as np
import pandas as pd
from audit_lib import (extract_pxd004242, load_pxd004242_meta, longitudinal, draw_one_run_per_sample,
                       score_matrix, pairs_from, auc, sens_at_spec, TIME_ORDER, WEEKS)

X = extract_pxd004242(); meta = longitudinal(load_pxd004242_meta())
rot = meta[meta.tech.astype(int) == ((meta.time.map(TIME_ORDER) - 1) % 4) + 1]


def prim(sel, features=None):
    S, f = score_matrix(X.loc[sel.index], sel.batch, features=features)
    P = pairs_from(S, sel, ("donor", "batch", "group"))
    return P[P.same_group & ~P.same_batch], f


# ---- claim 6
rows = []
for seed in range(1000, 1020):
    feats = sorted(np.random.default_rng(seed).choice(X.columns.to_numpy(), 273, replace=False).tolist())
    P, _ = prim(rot, feats)
    r = {"seed": seed, "rot_auc": auc(P.same_donor, P.score), "rot_sens99": sens_at_spec(P.same_donor, P.score)}
    sel = draw_one_run_per_sample(meta, np.random.default_rng(seed))
    P, _ = prim(sel, feats)
    r.update(rand_auc=auc(P.same_donor, P.score), rand_sens99=sens_at_spec(P.same_donor, P.score))
    cov = X.loc[meta.index, feats].where(X.loc[meta.index, feats] > 0).notna().mean()
    r["median_coverage"] = cov.median()
    rows.append(r)
R = pd.DataFrame(rows)
print("random-273 panels (20):")
print(R.drop(columns="seed").describe().loc[["50%", "min", "max"]].round(4).to_string())
cov_all = X.loc[meta.index].where(X.loc[meta.index] > 0).notna().mean()
print(f"features with >=99% coverage: {(cov_all >= .99).sum()} of {len(cov_all)}; median coverage all={cov_all.median():.3f}")

# ---- claim 5: per-gap AUC, 50 random assignments
recs = []
for d in range(50):
    sel = draw_one_run_per_sample(meta, np.random.default_rng(20260919 + d))
    P, _ = prim(sel)
    w = sel.time.map(WEEKS).to_numpy()
    P = P.assign(gap=np.abs(w[P.i] - w[P.j]), pre=(np.minimum(w[P.i], w[P.j]) < 0))
    neg = P.score[~P.same_donor].to_numpy()
    pos = P[P.same_donor]
    for g, b in pos.groupby("gap"):
        y = np.r_[np.ones(len(b)), np.zeros(len(neg))]
        recs.append({"draw": d, "gap": g, "n": len(b), "auc": auc(y, np.r_[b.score, neg])})
    for lab, b in (("pre_to_post", pos[pos.pre]), ("post_only", pos[~pos.pre])):
        y = np.r_[np.ones(len(b)), np.zeros(len(neg))]
        recs.append({"draw": d, "gap": lab, "n": len(b), "auc": auc(y, np.r_[b.score, neg])})
T = pd.DataFrame(recs)
G = T.groupby(T.gap.astype(str)).agg(n=("n", "median"), auc=("auc", "median"),
                                      lo=("auc", lambda v: v.quantile(.025)), hi=("auc", lambda v: v.quantile(.975)))
print("\nper-gap AUC (median over 50 assignments):"); print(G.round(4).to_string())
