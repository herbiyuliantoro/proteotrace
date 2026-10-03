"""Attack 1: optimism of tuning --margin on the cohorts it is reported on.

Nested design: split donors into halves A/B, run the WHOLE tool on each half as
its own cohort, tune the margin on A with planted swaps, report it on B, and
compare with B's own in-sample optimum. Plus cross-cohort transfer, plus the
rival-adjust check.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from a3_tool import load, score_once, verdicts, uniform_planter, repo_planter, account

GRID = np.round(np.arange(0.0, 0.501, 0.025), 3)


def sweep(X, S, n_rep, seed, planter=uniform_planter, n_pairs=12):
    scores = score_once(X, S)
    rows = []
    for r in range(n_rep):
        rng = np.random.default_rng(seed + r)
        t, sw = planter(S, n_pairs, rng)
        for m in GRID:
            a = account(verdicts(scores, t, m), S, t, sw)
            a.update(margin=m, rep=r)
            rows.append(a)
    D = pd.DataFrame(rows).groupby("margin").mean(numeric_only=True)
    D["J"] = D.sens - D.fa_repo_def
    return D


def pick(D):
    return float(D.J.idxmax())


rng = np.random.default_rng(7)
for name in ("PXD001064", "PXD004242"):
    X, S = load(name)
    if name == "PXD001064":
        unit = S.donor.str.extract(r"(Twin pair-\d+)")[0]
    else:
        unit = S.donor
    units = unit.unique()
    res = []
    for k in range(8):
        u = units.copy(); rng.shuffle(u)
        halves = [set(u[: len(u) // 2]), set(u[len(u) // 2:])]
        sub = [S[unit.isin(h)].reset_index(drop=True) for h in halves]
        D = [sweep(X.loc[s.run], s, 12, 900 + 31 * k + 7 * h, n_pairs=6) for h, s in enumerate(sub)]
        for a, b in ((0, 1), (1, 0)):
            m_a, m_b = pick(D[a]), pick(D[b])
            res.append({"split": k, "m_tuned": m_a, "m_oracle": m_b,
                        "J_tuned_on_A_seen_on_A": D[a].J[m_a],
                        "J_oracle_B": D[b].J[m_b], "J_tuned_A_applied_B": D[b].J[m_a],
                        "J_default_0.20_on_B": D[b].J[0.2],
                        "sens_B_at_tuned": D[b].sens[m_a], "fa_B_at_tuned": D[b].fa_repo_def[m_a],
                        "sens_B_at_0.20": D[b].sens[0.2], "fa_B_at_0.20": D[b].fa_repo_def[0.2]})
    R = pd.DataFrame(res)
    R["optimism_J"] = R.J_tuned_on_A_seen_on_A - R.J_tuned_A_applied_B
    print(f"\n#### {name}: nested donor-half tuning (8 splits x 2 directions, 12 swap repeats each)")
    print(R.drop(columns="split").describe().loc[["mean", "50%", "min", "max"]].round(3).to_string())

# cross-cohort: tune on one full cohort, apply to the other
full = {n: sweep(*load(n), 30, 4000) for n in ("PXD001064", "PXD004242")}
for a, b in (("PXD001064", "PXD004242"), ("PXD004242", "PXD001064")):
    m = pick(full[a])
    print(f"\ntune on {a} -> m*={m}; on {b}: sens={full[b].sens[m]:.3f} fa={full[b].fa_repo_def[m]:.3f} "
          f"J={full[b].J[m]:.3f} | {b}'s own optimum m={pick(full[b])} J={full[b].J.max():.3f}")
for n, D in full.items():
    print(f"\n{n} full-cohort sweep (uniform swaps, 30 repeats):")
    print(D[["sens", "fa_repo_def", "fa_uninvolved", "fa_collateral", "J"]].round(3).to_string())

# repo planter, same grid: where is "the best point"?
for n in ("PXD001064", "PXD004242"):
    D = sweep(*load(n), 5, 100, planter=repo_planter)
    print(f"\n{n} REPO planter seeds 100-104, J-optimal margin = {pick(D)}")
    print(D.loc[[0.05, 0.1, 0.15, 0.2, 0.25, 0.3], ["sens", "fa_repo_def", "J"]].round(3).to_string())

# rival-adjust (Bonferroni-style) vs margin at matched false-alarm rate
print("\n#### rival-adjust vs margin (uniform swaps, 20 repeats)")
for n in ("PXD001064", "PXD004242"):
    X, S = load(n)
    scores = score_once(X, S)
    rows = []
    for r in range(20):
        t, sw = uniform_planter(S, 12, np.random.default_rng(7000 + r))
        for adj in (1, 10, 100, 1000):
            for m in (0.0, 0.05, 0.1, 0.2, 0.3):
                a = account(verdicts(scores, t, m, rival_adjust=adj), S, t, sw)
                a.update(adj=adj, margin=m)
                rows.append(a)
    D = pd.DataFrame(rows).groupby(["adj", "margin"])[["sens", "fa_repo_def", "fa_uninvolved"]].mean()
    print(n); print(D.round(3).to_string())
