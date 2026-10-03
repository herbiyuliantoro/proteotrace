"""Claim 4: donor-cluster bootstrap CIs, and threshold portability across donors."""
import numpy as np
import pandas as pd
from audit_lib import (load_pxd001064, extract_pxd004242, load_pxd004242_meta, longitudinal,
                       draw_one_run_per_sample, score_matrix, pairs_from, auc, sens_at_spec)

X1, m1 = load_pxd001064()
S1, _ = score_matrix(X1, m1.batch)
P1 = pairs_from(S1, m1, ("donor", "batch", "pair_id"))
P1 = P1[~P1.same_batch]
u1 = m1.pair_id.to_numpy()

X2 = extract_pxd004242(); meta = longitudinal(load_pxd004242_meta())
sel = draw_one_run_per_sample(meta, np.random.default_rng(20261002))
S2, _ = score_matrix(X2.loc[sel.index], sel.batch)
P2 = pairs_from(S2, sel, ("donor", "batch", "group"))
P2 = P2[P2.same_group & ~P2.same_batch]
u2 = sel.donor.to_numpy()

for name, P, u in (("PXD001064", P1, u1), ("PXD004242 (one random draw)", P2, u2)):
    units = np.unique(u)
    rng = np.random.default_rng(0)
    a, s = [], []
    ui, uj = u[P.i], u[P.j]
    for b in range(300):
        w = pd.Series(rng.multinomial(len(units), np.ones(len(units)) / len(units)), index=units)
        wt = w.loc[ui].to_numpy() * w.loc[uj].to_numpy()       # pair weight = product of cluster counts
        keep = wt > 0
        idx = np.repeat(np.arange(len(P))[keep], wt[keep])
        Q = P.iloc[idx]
        a.append(auc(Q.same_donor, Q.score)); s.append(sens_at_spec(Q.same_donor, Q.score))
    print(f"{name}: AUC={auc(P.same_donor, P.score):.4f} cluster-bootstrap 95% [{np.quantile(a,.025):.4f},{np.quantile(a,.975):.4f}]"
          f" | sens99={sens_at_spec(P.same_donor, P.score):.3f} [{np.quantile(s,.025):.3f},{np.quantile(s,.975):.3f}]")
    # threshold portability: 99% quantile from half the donors' negatives, applied to the other half
    sp, se = [], []
    for k in range(200):
        r = np.random.default_rng(k)
        test = set(r.choice(units, len(units) // 2, replace=False))
        ti, tj = np.isin(ui, list(test)), np.isin(uj, list(test))
        tr, te = P[~ti & ~tj], P[ti & tj]
        thr = np.quantile(tr.score[~tr.same_donor], .99)
        sp.append(np.mean(te.score[~te.same_donor] < thr)); se.append(np.mean(te.score[te.same_donor] >= thr))
    print(f"   threshold from other donors -> realised specificity median {np.median(sp):.4f} "
          f"[{np.quantile(sp,.025):.4f},{np.quantile(sp,.975):.4f}], sensitivity median {np.median(se):.3f}")
