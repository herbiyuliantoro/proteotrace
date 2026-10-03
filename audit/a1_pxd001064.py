"""Claims 1-4 on PXD001064, re-derived independently."""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from audit_lib import load_pxd001064, score_matrix, pairs_from, auc, sens_at_spec

for dedup in ("first", "last"):
    X, meta = load_pxd001064(dedup)
    print(f"\n=== dedup={dedup}: runs={len(meta)} donors={meta.donor.nunique()} batches={meta.batch.nunique()}")
    print("batch sizes", meta.batch.value_counts().sort_index().to_dict())
    S_raw, feats = score_matrix(X, meta.batch, correct=False, scale=False)
    S_cor, _ = score_matrix(X, meta.batch)
    print("features", len(feats), "contaminant features:", [f for f in feats if "CONT_" in f])
    for name, S in (("raw", S_raw), ("corrected", S_cor)):
        p = pairs_from(S, meta)
        cb = p[~p.same_batch]
        print(f"{name:9s} cross-batch n={len(cb)} pos={cb.same_donor.sum()} "
              f"AUC={auc(cb.same_donor, cb.score):.4f} "
              f"sens95={sens_at_spec(cb.same_donor, cb.score, .95):.4f} "
              f"sens99={sens_at_spec(cb.same_donor, cb.score, .99):.4f}")
    # raw *with* z-scaling but no batch correction (isolates the correction step)
    S_rz, _ = score_matrix(X, meta.batch, correct=False, scale=True)
    p = pairs_from(S_rz, meta); cb = p[~p.same_batch]
    print(f"raw+zscale cross-batch AUC={auc(cb.same_donor, cb.score):.4f}")

    # ladder
    p = pairs_from(S_cor, meta, cols=("donor", "batch", "pair_id"))
    zyg = meta["zygosity"].to_numpy()
    rel = np.where(p.same_donor, "same", np.where(p.same_pair_id,
                   np.where(zyg[p.i] == "MZ", "MZ", "DZ"), "unrelated"))
    # labelling sanity: both members of a twin pair must share zygosity
    zz = meta.groupby("pair_id")["zygosity"].nunique()
    print("pairs with inconsistent zygosity:", list(zz[zz > 1].index))
    print("twin pairs MZ/DZ:", meta.groupby("pair_id")["zygosity"].first().value_counts().to_dict())
    p["rel"] = rel
    cb = p[~p.same_batch]
    for r in ("same", "MZ", "DZ", "unrelated"):
        v = cb.loc[cb.rel == r, "score"]
        print(f"  {r:9s} n={len(v):5d} median={v.median():.4f}")
    for a, b in (("same", "MZ"), ("MZ", "DZ"), ("DZ", "unrelated"), ("same", "unrelated")):
        x = cb.loc[cb.rel == a, "score"]; y = cb.loc[cb.rel == b, "score"]
        u = mannwhitneyu(x, y).statistic
        print(f"  P({a}>{b})={u/(len(x)*len(y)):.4f}")
    # estimand without twins in the negative class (comparable to PXD004242, no twins)
    nt = cb[cb.rel.isin(["same", "unrelated"])]
    print(f"AUC same vs unrelated only = {auc(nt.rel == 'same', nt.score):.4f}, "
          f"sens99 = {sens_at_spec(nt.rel == 'same', nt.score):.4f}")
    # donor-level dependence: number of positives per donor
    print("positives per donor:", cb[cb.same_donor].groupby(meta.donor.to_numpy()[cb[cb.same_donor].i]).size().value_counts().to_dict())
