"""Fig 6B: closed-set, in-group identification vs the open problem."""
import numpy as np
import pandas as pd
from audit_lib import extract_pxd004242, load_pxd004242_meta, longitudinal, draw_one_run_per_sample, score_matrix

X = extract_pxd004242()
meta = longitudinal(load_pxd004242_meta())
res = []
for d in range(5):
    sel = draw_one_run_per_sample(meta, np.random.default_rng(300 + d))
    S, _ = score_matrix(X.loc[sel.index], sel.batch)
    donor = sel.donor.to_numpy(); group = sel.group.to_numpy(); plate = sel.batch.to_numpy()
    n = len(sel)

    def top(q, excluded, in_group):
        best, bd = -np.inf, None
        for dn in np.unique(donor):
            m = (donor == dn)
            m[list(excluded)] = False
            m[q] = False
            if in_group:
                m &= group == group[q]
            if m.sum() == 0:
                continue
            v = np.nanmedian(S[q, m])
            if v > best:
                best, bd = v, dn
        return bd

    # clean runs: how often does the top donor disagree with a CORRECT label?
    for ig in (True, False):
        dis = np.mean([top(q, set(), ig) != donor[q] for q in range(n)])
        res.append({"draw": d, "candidates": "in-group" if ig else "all 43", "what": "clean disagreement", "rate": dis})
    rng = np.random.default_rng(900 + d)
    swaps = []
    while len(swaps) < 200:
        i, j = rng.choice(n, 2, replace=False)
        if plate[i] == plate[j] and donor[i] != donor[j]:
            swaps.append((i, j))
    for ig in (True, False):
        for excl in (True, False):
            rec = dis = 0
            for i, j in swaps:
                for q, assigned in ((i, donor[j]), (j, donor[i])):
                    t = top(q, {i, j} if excl else set(), ig)
                    dis += t != assigned
                    rec += t == donor[q]
            res.append({"draw": d, "candidates": "in-group" if ig else "all 43",
                        "what": f"swap recovery, partner {'excluded' if excl else 'kept'}", "rate": rec / 400})
            res.append({"draw": d, "candidates": "in-group" if ig else "all 43",
                        "what": f"swap disagreement, partner {'excluded' if excl else 'kept'}", "rate": dis / 400})
R = pd.DataFrame(res)
print(R.groupby(["what", "candidates"]).rate.mean().round(4).to_string())
