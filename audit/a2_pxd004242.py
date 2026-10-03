"""Claims 1,2,4,5,6 on PXD004242, re-derived from the MaxQuant zip + SDRF."""
import sys
import numpy as np
import pandas as pd
from audit_lib import (ROOT, extract_pxd004242, load_pxd004242_meta, longitudinal,
                       score_matrix, pairs_from, auc, sens_at_spec, TIME_ORDER, WEEKS)

X = extract_pxd004242()
tracked = pd.read_csv(f"{ROOT}/data/PXD004242/PXD004242_protein_matrix.tsv", sep="\t", index_col=0)
same_shape = X.shape == tracked.shape
Xa = X.loc[tracked.index, tracked.columns]
print("extracted", X.shape, "tracked", tracked.shape,
      "identical:", same_shape and np.allclose(Xa.fillna(-1).to_numpy(), tracked.fillna(-1).to_numpy()))

meta_all = load_pxd004242_meta()
print("label donor == SDRF individual:", (meta_all.label_donor == meta_all.donor).mean(),
      "| label tidx == time order:",
      (meta_all.label_tidx.astype(int) == meta_all.time.map(TIME_ORDER)).mean())
meta = longitudinal(meta_all)
print("longitudinal donors", meta.donor.nunique(), "samples", meta["sample"].nunique(), "runs", len(meta))
print("donors per plate group:", meta.groupby("group").donor.nunique().to_dict())
print("donor in >1 group:", (meta.groupby("donor").group.nunique() > 1).sum())

# Cross-check against the repo's derived table (structure only)
bio = pd.read_csv(f"{ROOT}/data/PXD004242/PXD004242_biological_samples.tsv", sep="\t", dtype=str)
bio = bio[bio.cohort == "longitudinal"]
rep_runs = set(";".join(bio.maxquant_labels).split(";"))
print("repo run set == SDRF-derived run set:", rep_runs == set(meta.index))
# technical replicate index in repo table vs SDRF
tmap = {}
for _, r in bio.iterrows():
    for lab, t in zip(r.maxquant_labels.split(";"), r.technical_replicates.split(",")):
        tmap[lab] = t
print("tech index repo == SDRF:", np.mean([tmap[r] == meta.loc[r, "tech"] for r in meta.index]))


def metrics(sel, features=None, tag=""):
    Xs = X.loc[sel.index]
    out = {}
    for name, kw in (("raw", dict(correct=False, scale=False)),
                     ("raw_z", dict(correct=False, scale=True)),
                     ("corr", dict())):
        S, feats = score_matrix(Xs, sel.batch, features=features, **kw)
        p = pairs_from(S, sel, cols=("donor", "batch", "group"))
        prim = p[p.same_group & ~p.same_batch]
        allx = p[~p.same_batch]
        out[name] = auc(prim.same_donor, prim.score)
        if name == "corr":
            out["sens95"] = sens_at_spec(prim.same_donor, prim.score, .95)
            out["sens99"] = sens_at_spec(prim.same_donor, prim.score, .99)
            out["all_xplate_auc"] = auc(allx.same_donor, allx.score)
            out["all_xplate_sens99"] = sens_at_spec(allx.same_donor, allx.score, .99)
            out["nfeat"] = len(feats)
            out["n_pos"] = int(prim.same_donor.sum())
            out["n_neg_within"] = int((~prim.same_donor).sum())
            out["neg_within_median"] = float(prim.score[~prim.same_donor].median())
            xg = p[~p.same_group]
            out["neg_crossgroup_median"] = float(xg.score.median())
            out["_pairs"] = (p, sel)
    return out


# --- rotation (deterministic)
rot = meta[meta.tech.astype(int) == ((meta.time.map(TIME_ORDER) - 1) % 4) + 1]
# participant 44 has 6 timepoints, still 4 runs per sample
assert rot["sample"].is_unique and len(rot) == meta["sample"].nunique()
r = metrics(rot)
print(f"\nROTATION raw={r['raw']:.4f} raw+z={r['raw_z']:.4f} corrected={r['corr']:.4f} "
      f"sens95={r['sens95']:.4f} sens99={r['sens99']:.4f} feats={r['nfeat']} pos={r['n_pos']} "
      f"neg={r['n_neg_within']} | all cross-plate AUC={r['all_xplate_auc']:.4f} sens99={r['all_xplate_sens99']:.4f}"
      f" | neg median within={r['neg_within_median']:.4f} cross-group={r['neg_crossgroup_median']:.4f}")

# --- random assignments, replicating the repo's draw procedure from its documented rule:
# samples in sorted sample_id order, rng.integers over the runs listed for that sample.
cand = {}
for _, row in bio.iterrows():
    cand[row.sample_id] = row.maxquant_labels.split(";")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
rows = []
top273 = None
cov = X.loc[meta.index].where(X.loc[meta.index] > 0).notna().mean(axis=0)
top273 = cov.sort_values(ascending=False, kind="mergesort").index[:273].tolist()
for k in range(N):
    rng = np.random.default_rng(20260919 + k)
    pick = [cand[s][int(rng.integers(0, len(cand[s])))] for s in sorted(cand)]
    sel = meta.loc[pick]
    m = metrics(sel)
    t = metrics(sel, features=top273)
    rows.append({"draw": k, "raw": m["raw"], "raw_z": m["raw_z"], "corr": m["corr"],
                 "sens95": m["sens95"], "sens99": m["sens99"], "nfeat": m["nfeat"],
                 "allx_auc": m["all_xplate_auc"], "allx_sens99": m["all_xplate_sens99"],
                 "t273_corr": t["corr"], "t273_sens99": t["sens99"], "t273_raw": t["raw"],
                 "n_pos": m["n_pos"]})
    if k % 50 == 0:
        print("draw", k, flush=True)
D = pd.DataFrame(rows)
D.to_csv("pxd004242_draws.tsv", sep="\t", index=False)
q = lambda c: f"{D[c].median():.4f} [{D[c].quantile(.025):.4f},{D[c].quantile(.975):.4f}]"
for c in D.columns[1:]:
    print(f"{c:12s} median [2.5,97.5] = {q(c)}")
