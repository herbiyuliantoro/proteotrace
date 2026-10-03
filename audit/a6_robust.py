"""Code-level checks of the shipped tool and known-answer cohorts it was not tested on."""
import sys, subprocess, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))
import numpy as np
import pandas as pd
from a3_tool import load, score_once, verdicts
from proteotrace.core import normalise, similarity

X, S = load("PXD001064")
ann = S.set_index("run")
scores = score_once(X, S)
base = verdicts(scores, S, 0.20)
print("clean PXD001064 verdicts:", base.verdict.value_counts().to_dict())

# 1. why no_reference?
nr = base[base.verdict == "no_reference"]
runs_per_donor = S.groupby("donor").run.size()
batches_per_donor = S.groupby("donor").batch.nunique()
d = ann.loc[nr.index, "donor"]
print(f"no_reference runs={len(nr)}: donor has 1 run: {(runs_per_donor[d] == 1).sum()}, "
      f"donor has >=2 runs all in ONE batch: {((runs_per_donor[d] >= 2) & (batches_per_donor[d] == 1)).sum()}")

# 2. row-order invariance
perm = np.random.default_rng(0).permutation(len(S))
S2 = S.iloc[perm].reset_index(drop=True)
X2 = X.iloc[np.random.default_rng(1).permutation(len(X))]
v2 = verdicts(score_once(X2.loc[S2.run], S2), S2, 0.20)
print("row-order invariant:", (v2.loc[base.index, "verdict"] == base.verdict).all())

# 3. a run that failed (95% of values missing): what does it get called?
Xf = X.copy()
victim = base.index[base.verdict == "ok"][0]
cols = Xf.columns[np.random.default_rng(2).random(Xf.shape[1]) < 0.95]
Xf.loc[victim, cols] = np.nan
vf = verdicts(score_once(Xf, S), S, 0.20)
print(f"failed run {victim}: verdict={vf.loc[victim, 'verdict']} n_ref={vf.loc[victim, 'n_reference_runs']} "
      f"features now={int(normalise(Xf, ann['batch']).n_features)} (was 273)")

# 4. tiny batches: one run alone in a batch; two runs alone in a batch
Sb = S.copy()
lone = base.index[base.verdict == "ok"][1]
Sb.loc[Sb.run == lone, "batch"] = "LONE"
duo = list(base.index[base.verdict == "ok"][2:4])
Sb.loc[Sb.run.isin(duo), "batch"] = "DUO"
sc = score_once(X, Sb)
vb = verdicts(sc, Sb, 0.20)
print(f"run alone in its batch: verdict={vb.loc[lone, 'verdict']} score_to_own={vb.loc[lone, 'score_to_claimed_donor']:.3f}"
      f" (was {base.loc[lone, 'score_to_claimed_donor']:.3f}); max |r| to anything={np.nanmax(np.abs(sc.loc[lone].drop(lone))):.3f}")
print(f"two runs alone in a batch: r between them={sc.loc[duo[0], duo[1]]:.3f}; verdicts={vb.loc[duo, 'verdict'].tolist()}")

# 5. 'review' = weak run? compare missingness of review vs ok runs (clean cohort)
miss = X.where(X > 0).isna().mean(axis=1)
print("missing fraction (all 712 groups) median: ok={:.3f} review={:.3f}".format(
    miss[base.index[base.verdict == "ok"]].median(), miss[base.index[base.verdict == "review"]].median()))
tic = np.log10(X.where(X > 0).sum(axis=1))
print("log10 total intensity median: ok={:.3f} review={:.3f}".format(
    tic[base.index[base.verdict == "ok"]].median(), tic[base.index[base.verdict == "review"]].median()))

# 6. known-answer cohorts the author did not test
def summary(v, mask=None):
    v = v if mask is None else v[mask]
    return v.verdict.value_counts().to_dict()

for name in ("PXD001064", "PXD004242"):
    X, S = load(name); sc = score_once(X, S)
    rng = np.random.default_rng(11)
    # (a) every label scrambled
    Sa = S.copy(); Sa["donor"] = rng.permutation(Sa.donor.to_numpy())
    print(f"\n{name} (a) ALL labels permuted:", summary(verdicts(sc, Sa, 0.20)))
    # (b) off-by-one: in one batch, every run gets its neighbour's label (a plate shift)
    Sb = S.copy()
    bname = Sb.batch.value_counts().index[0]
    idx = Sb.index[Sb.batch == bname]
    Sb.loc[idx, "donor"] = np.roll(Sb.loc[idx, "donor"].to_numpy(), 1)
    shifted = Sb.run[idx][(Sb.loc[idx, "donor"].to_numpy() != S.loc[idx, "donor"].to_numpy())]
    vb = verdicts(sc, Sb, 0.20)
    print(f"{name} (b) off-by-one shift in batch {bname} ({len(shifted)} runs relabelled):",
          summary(vb, vb.index.isin(shifted)), "| rest:", summary(vb, ~vb.index.isin(shifted)))
    # (c) a sample from someone NOT in the cohort: drop donor D, relabel one of D's runs as donor B,
    #     and drop B's run on that batch so B still has its other runs as reference.
    out = []
    for t in range(40):
        r = np.random.default_rng(100 + t)
        donors = S.donor.unique()
        D, B = r.choice(donors, 2, replace=False)
        intruder = S[S.donor == D].run.iloc[0]
        Sc = S[(S.donor != D) | (S.run == intruder)].copy()
        Sc.loc[Sc.run == intruder, "donor"] = B
        Sc = Sc.reset_index(drop=True)
        v = verdicts(score_once(X.loc[Sc.run], Sc), Sc, 0.20)
        out.append(v.loc[intruder, "verdict"])
    print(f"{name} (c) intruder from outside the cohort, 40 trials:", pd.Series(out).value_counts().to_dict())
    # (d) degraded data: multiplicative noise sd 0.5 log2 on every value
    Xd = X * np.exp2(np.random.default_rng(3).normal(0, 0.5, X.shape))
    print(f"{name} (d) noisy clean cohort:", summary(verdicts(score_once(Xd, S), S, 0.20)))
    # (e) batch column omitted
    Se = S.copy(); Se["batch"] = "single_batch"
    ve = verdicts(score_once(X, Se), Se, 0.20)
    print(f"{name} (e) no batch annotation:", summary(ve))

# 7. the CLI end-to-end reproduces the tutorial
out = subprocess.run([sys.executable, "-m", "proteotrace.cli", "check",
                      "--matrix", "data/PXD001064_protein_matrix.tsv",
                      "--samples", "examples/PXD001064_samples.tsv", "--out", sys.argv[1]],
                     cwd=str(__import__("pathlib").Path(__file__).resolve().parent.parent), capture_output=True, text=True)
print("\nCLI:", out.stdout[:400])
