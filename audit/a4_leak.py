"""Attack 2 (held-out fitting), step attribution, attack 5 (residual plate signal)."""
import numpy as np
import pandas as pd
from audit_lib import (load_pxd001064, extract_pxd004242, load_pxd004242_meta, longitudinal,
                       score_matrix, pairs_from, auc, sens_at_spec, TIME_ORDER)

# ------------------------------------------------------------------ cohorts
X1, m1 = load_pxd001064()
m1["unit"] = m1["pair_id"]                     # split by twin pair: twins never straddle
X2 = extract_pxd004242()
meta2 = longitudinal(load_pxd004242_meta())
rot = meta2[meta2.tech.astype(int) == ((meta2.time.map(TIME_ORDER) - 1) % 4) + 1].copy()
rot["unit"] = rot["donor"]
X2r = X2.loc[rot.index]


def eval_mask(p, cohort):
    if cohort == 1:
        return ~p.same_batch
    return p.same_group & ~p.same_batch


def cols(cohort):
    return ("donor", "batch") if cohort == 1 else ("donor", "batch", "group")


print("=== step attribution (AUC on the primary estimand)")
for cohort, X, m in ((1, X1, m1), (2, X2r, rot)):
    out = {}
    for name, kw in (("raw", dict(correct=False, scale=False)),
                     ("zscale only", dict(correct=False, scale=True)),
                     ("batch only", dict(correct=True, scale=False)),
                     ("batch + z", dict())):
        S, _ = score_matrix(X, m.batch, **kw)
        p = pairs_from(S, m, cols(cohort)); e = p[eval_mask(p, cohort)]
        out[name] = f"{auc(e.same_donor, e.score):.4f}/{sens_at_spec(e.same_donor, e.score):.3f}"
    print(f"  {'PXD001064' if cohort == 1 else 'PXD004242 rotation'}: AUC/sens99", out)

if True:
    # contaminant-tagged features in PXD001064
    S, feats = score_matrix(X1, m1.batch)
    pure = [f for f in feats if all("CONT_" in t for t in f.split("/")[1:])]
    anyc = [f for f in feats if "CONT_" in f]
    for tag, drop in (("drop pure contaminants", pure), ("drop any CONT-tagged group", anyc)):
        keep = [f for f in feats if f not in drop]
        S, _ = score_matrix(X1, m1.batch, features=keep)
        p = pairs_from(S, m1); e = p[~p.same_batch]
        print(f"  PXD001064 {tag} ({len(drop)} removed, {len(keep)} left): "
              f"AUC={auc(e.same_donor, e.score):.4f} sens99={sens_at_spec(e.same_donor, e.score):.4f}")

print("\n=== attack 2: fit on training donors, score held-out donors")
rng = np.random.default_rng(1)
for cohort, X, m in ((1, X1, m1), (2, X2r, rot)):
    rows = []
    S_full, _ = score_matrix(X, m.batch)
    P_full = pairs_from(S_full, m, cols(cohort))
    units = m["unit"].unique()
    for k in range(50):
        if cohort == 2:   # stratify the split by plate group
            test_units = set()
            for g, b in m.groupby("group"):
                u = b["unit"].unique(); rng.shuffle(u); test_units |= set(u[: len(u) // 2])
        else:
            u = units.copy(); rng.shuffle(u); test_units = set(u[: len(u) // 2])
        test = m["unit"].isin(test_units).to_numpy()
        S_ho, _ = score_matrix(X, m.batch, fit=~test)
        P_ho = pairs_from(S_ho, m, cols(cohort))
        def sub(P):
            tt = test[P.i] & test[P.j]
            tr = ~test[P.i] & ~test[P.j]
            em = eval_mask(P, cohort)
            return P[tt & em], P[tr & em]
        te_full, _ = sub(P_full)
        te_ho, tr_ho = sub(P_ho)
        # threshold set on training pairs, applied to held-out pairs
        thr = np.quantile(tr_ho.score[~tr_ho.same_donor], .99)
        rows.append({
            "auc_transductive": auc(te_full.same_donor, te_full.score),
            "auc_heldout": auc(te_ho.same_donor, te_ho.score),
            "sens99_insample_thr": sens_at_spec(te_ho.same_donor, te_ho.score),
            "sens99_train_thr": float(np.mean(te_ho.score[te_ho.same_donor] >= thr)),
            "spec_train_thr": float(np.mean(te_ho.score[~te_ho.same_donor] < thr)),
            "n_pos_test": int(te_ho.same_donor.sum()),
        })
    D = pd.DataFrame(rows)
    D["drop"] = D.auc_transductive - D.auc_heldout
    print(f"  {'PXD001064' if cohort == 1 else 'PXD004242 rotation'} (50 donor splits): ")
    print("   ", D.describe().loc[["mean", "50%", "min", "max"]].round(4).to_string().replace("\n", "\n    "))

print("\n=== attack 5: residual plate / date structure in PXD004242 (rotation)")
S, _ = score_matrix(X2r, rot.batch)
P = pairs_from(S, rot, ("donor", "batch", "group", "date"))
prim = P[P.same_group & ~P.same_batch].copy()
plates = rot.batch.to_numpy()
prim["pp"] = [tuple(sorted(t)) for t in zip(plates[prim.i], plates[prim.j])]
# stratified AUC: positives compared only with negatives from the same plate pair
num = den = 0.0
for _, g in prim.groupby("pp"):
    if g.same_donor.any() and (~g.same_donor).any():
        w = g.same_donor.sum() * (~g.same_donor).sum()
        num += auc(g.same_donor, g.score) * w; den += w
print(f"  pooled primary AUC={auc(prim.same_donor, prim.score):.4f}  plate-pair-stratified AUC={num/den:.4f}")
neg = prim[~prim.same_donor]
print(f"  within-group cross-plate negatives: same acquisition date n={neg.same_date.sum()} "
      f"median={neg.score[neg.same_date].median():.4f} | different date median={neg.score[~neg.same_date].median():.4f}")
print(f"  share-date rate: positives={prim.same_date[prim.same_donor].mean():.3f} negatives={neg.same_date.mean():.3f}")
# can the corrected score still tell which plate group a pair comes from?
dd = P[~P.same_donor & ~P.same_batch]
print(f"  different-donor cross-plate: within-group median={dd.score[dd.same_group].median():.4f} "
      f"cross-group median={dd.score[~dd.same_group].median():.4f}; "
      f"AUC(within vs cross group)={auc(dd.same_group, dd.score):.4f}")
# same plate, different donor vs different plate same group
sp = P[~P.same_donor & P.same_batch]
print(f"  different-donor SAME plate median={sp.score.median():.4f} (n={len(sp)})")
