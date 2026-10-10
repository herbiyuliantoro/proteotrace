# ProteoTrace: analysis and software audit

Audit of branch `audit` at `323818e` (= `main` at `5f5d5b4` plus the original
audit brief). Audit date: 2026-10-02. See the [audit scope](AUDIT_SCOPE.md).

All numbers below come from scripts in [`audit/`](../audit). The scoring code
there (`audit/audit_lib.py`) is a fresh implementation written from the
method description. It reads only the raw inputs: the PXD001064 matrix and
the PRIDE sample map, and for PXD004242 the MaxQuant zip and community SDRF.
It never imports the repo's analysis code or reads `data/*.json`. The tool
checks (`a3`, `a5`, `a6`) do call the shipped `proteotrace` package, because
claim 7 is a claim about that package. They plant swaps and keep the
accounting independently. Re-run everything with `cd audit && sh run_all.sh`;
saved outputs are in `audit/outputs/`.

Those outputs were captured against `proteotrace` 0.1.0, i.e. *before* the
fixes in section 5. Re-running `a3`, `a5` and `a6` against 0.2.0 shows the
fixed behaviour. One example: the failed run in F10 now gets `unscorable`.

---

## Audit summary

The central discrimination result was reproduced using a separate
implementation and the raw inputs:

- Both AUCs reproduce to four decimals.
- Fitting the preprocessing on held-out donors costs only 0.005–0.015 AUC.
- No residual plate-pair leakage was detected in the checks performed.

The audit identified two problems with the interpretation and validation of
the version examined:

1. **"Batch correction is the method" is not what the data show.** The
   "raw" baseline (0.7827 / 0.8463) leaves out *two* steps: the batch-median
   correction and the per-feature z-scaling. Z-scaling alone takes
   PXD001064 to 0.9313 and PXD004242 to 0.9618. Batch correction then adds
   +0.038 and +0.012. The step from 0.78 to 0.97 belongs to the two steps
   together, and z-scaling is the larger share of it.
2. **Claim 7 has three problems:**
   - **One batch.** The swaps are planted almost entirely in the
     alphabetically first batch, in all five repeats.
   - **"False alarms" mostly counts collateral.** In PXD001064 they are
     mostly the *partner runs* of the swapped donors. Those runs are
     necessarily flagged when a donor has two visits.
   - **`review` is left out.** The table ignores the 17–23% of clean runs
     that land in `review`.

   Planting swaps uniformly across batches gives **72–74% / 97–98%**
   sensitivity, against the reported 80.6% / 93.3%.

Nested margin tuning costs 0.01–0.02 in Youden J. A margin tuned on
either cohort lands within 0.008 J of the other cohort's own optimum.

---

## 1. Verdict per claim

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1 | AUC 0.9692 / 0.9733 | **Holds, with caveat** | Reproduced exactly: 0.9692, and 0.9733 [0.9633, 0.9804] over 200 draws (`a1`, `a2`). Donor-cluster bootstrap 95% CIs are [0.937, 0.986] and [0.949, 0.987] (`a8`). The two estimands differ (F5). |
| 2 | Batch correction is the method: 0.78→0.97, 0.85→0.97 | **Does not hold as worded** | Numbers reproduce, but the baseline also omits z-scaling. Z-scaling alone gives 0.9313 / 0.9618; batch correction then adds +0.038 / +0.012 (F1). |
| 3 | Ladder 0.381 > 0.293 > 0.156 > −0.007 | **Holds** | Reproduced exactly (`a1`). Zygosity is consistent within all 58 pairs; DZ and unrelated are disjoint by construction (different `Twin pair-NN`). The limit is understated (F9). |
| 4 | Sensitivity 0.677 / 0.697 at 99% specificity | **Holds, with caveat** | Reproduced: 0.6774 (63/93), and 0.6967 [0.648, 0.746] over 200 draws. Donor-cluster CIs are [0.42, 0.80] and [0.52, 0.80] (`a8`). The quoted assignment ranges understate uncertainty several-fold. A threshold taken from other donors' pairs gives realised specificity 0.9896 in both cohorts. |
| 5 | Survives intervention; gap-wise AUC 0.96–0.99 | **Holds, with caveat** | Re-derived over 50 draws: per-gap medians 0.956–0.988, pre→post 0.964, post-only 0.977 (`a9`). Each gap has about 32 positive pairs (about one per donor) and every per-gap range overlaps. "12% weight loss" is the source paper's figure, not measured here. |
| 6 | Top-273 0.976/0.715 vs random-273 0.953/0.526 | **Holds, with caveat** | Top-273 reproduced exactly (`a2`). Random-273 reproduced exactly under the rotation assignment (`a9`). Under the same random assignments as top-273, random-273 gives 0.9347/0.466, so the published contrast is *conservative*. The random panels have median coverage 0.25, so this compares dense with sparse features, not "quality vs size" in general. |
| 7 | 93.3% at 4.9%; 80.6% at 12.2% | **Does not hold as a description of the tool** | Reproduced with the repo's script. The planted set is one batch, the "false alarm" figure mixes in collateral, and `review` is omitted (F2). Representative numbers are in F2. |

Results by audit question:

| Audit question | Result |
|---|---|
| 1. Margin tuned on the reported cohorts | Real but small: about 0.01–0.02 in J (F3). |
| 2. Nothing held out | Unsupervised steps leak little: −0.005 / −0.015 AUC. The threshold transfers across donors. Batch medians must come from the cohort itself (F7). |
| 3. Different estimands | Confirmed. Matched estimand: 0.9707 vs 0.9674 (F5). |
| 4. Permutation vs physical swap | Faithful for reciprocal within-plate swaps. Optimistic for intruders, where the tool names the wrong donor (F6). |
| 5. Nested plate group | No plate-pair leak. Plate-group structure does remain in the score and flatters the within-group estimand (F5). |

---

## 2. Findings, ordered by how much they change the paper

### F1. The raw→corrected contrast credits batch correction with what z-scaling does

**What is wrong.** The baseline in every analysis script is log2,
run-median centring and Pearson. The "corrected" score adds *two* steps:
per-batch median subtraction, then per-feature z-scaling. All the scripts do
this (`analysis_cross_batch.py`, `residual_scores()` in
`analysis_pxd004242_assignment_sensitivity.py`, `proteotrace/core.py`). The
README, Figure 2's panel title ("After per-batch median subtraction") and
`ANALYSIS_LOG.md` credit the whole gain to batch correction.

**Demonstration.** `audit/a4_leak.py` runs a 2×2 ablation (AUC /
sensitivity at 99% specificity):

```text
PXD001064:          raw 0.7827/0.118  zscale only 0.9313/0.409  batch only 0.8815/0.215  batch + z 0.9692/0.677
PXD004242 rotation: raw 0.8425/0.258  zscale only 0.9704/0.636  batch only 0.9656/0.636  batch + z 0.9839/0.776
```

Over the 200 random draws (`audit/a2_pxd004242.py`): raw 0.8463, z-scale
only **0.9618** [0.9512, 0.9708], both 0.9733.

**How much it moves.** For AUC, batch correction's own contribution is
+0.038 (PXD001064) and +0.012 to +0.014 (PXD004242), not +0.19 and +0.13.
At 99% specificity the two steps interact strongly:

- PXD001064: z-scaling alone reaches 0.409; adding batch correction brings
  it to 0.677.
- PXD004242: 0.636 becomes 0.776.

So batch correction matters most exactly where the tool operates. That is a
defensible claim; "batch correction is the whole method" is not.

**Recommended change.** Put the 2×2 table in the paper (Table 1 or a
supplement). Retitle Figure 2 so it names both steps, or add the z-only
column. Rewrite the README's "The trick is..." paragraph to say that
per-protein standardisation plus within-batch centring together turn a weak
comparison into a usable one, with batch centring carrying the
high-specificity end.

### F2. Claim 7's swap set is one batch, its "false alarms" are mostly collateral, and `review` is unreported

**What is wrong.**

1. **One batch.** `validate_tool.py::plant_swaps` walks
   `groupby("batch")` in sorted order and stops after 12 pairs. Because
   batches hold 15–36 runs, every repeat swaps runs from the first batch.
   The rng only reshuffles *which* runs in that batch.
2. **Collateral counted as false alarms.** With two visits per donor, a
   swap between A and B leaves A's and B's *other* runs matched against a
   wrong reference. Those partner runs are then flagged by construction.
   The README calls the result "Flags clean runs".
3. **`review` left out.** Only `swap_suspected` is counted as a flag. On a
   clean cohort, 17–23% of runs come back `review`.

**Demonstration** (`audit/a3_tool.py`, the shipped `run_report`):

```text
PXD001064 repo planter seed 100: swapped runs by batch {'130104': 24}   (seeds 101-104 the same, one stray pair in 101)
PXD004242 repo planter seed 100: swapped runs by batch {'1B': 22, '1C': 2}  (all five seeds: 1B/1C only)
repo planter m=0.20 decomposed, PXD001064: sens 0.806, fa_repo_def 0.122,
      fa_collateral 0.837 (other runs of swapped donors), fa_uninvolved 0.034
CLEAN cohort, m=0.20: PXD001064 {'ok':126,'review':53,'no_reference':45,'swap_suspected':7}
                      PXD004242 {'ok':243,'review':52,'swap_suspected':5}
```

Uniform within-batch swaps over all batches. Twelve pairs per repeat, no
donor in two swaps, margin 0.20. Runs of 40 and 30 repeats with different
seeds (`a3`, `a5`):

| | sensitivity | "false alarm" (repo definition) | uninvolved runs flagged | collateral runs flagged |
|---|---|---|---|---|
| PXD001064 | 0.726 / 0.735 | 0.125 / 0.126 | 0.036 / 0.040 | 0.787 / 0.769 |
| PXD004242 | 0.967 / 0.982 | 0.034 / 0.032 | 0.017 / 0.015 | 0.049 / 0.048 |

**How much it moves.** PXD001064 sensitivity is overstated by about 7–8
points (80.6 → 73), and PXD004242 is *understated* by 3–5 points
(93.3 → 97–98). The 12.2% false-alarm figure becomes "about 4% of
uninvolved runs, plus about 80% of the swapped donors' other runs". The
second part is not a false alarm. It is the tool saying that one of four
runs is wrong without being able to say which.

**Recommended change.** Plant swaps uniformly across batches. Use 30 or
more repeats and report Monte-Carlo spread. Report three rates: planted
runs caught, uninvolved runs flagged, and other runs of affected donors
flagged. Add the clean-cohort verdict table (including `review`) to the
README and tutorial.

### F3. Margin-tuning optimism (audit question 1)

**Demonstration** (`audit/a5_margin.py`). Each cohort is split into donor
halves, by twin pair for PXD001064. The tool is run on each half as its own
cohort. The margin is tuned on half A by Youden J (sensitivity minus the
repo-definition false-alarm rate) over a 0–0.5 grid in 0.025 steps, then
reported on half B. 8 splits × 2 directions:

```text
                       J seen on A  J of A's margin on B  B's own optimum  optimism (mean / median)
PXD001064                   0.497          0.481              0.497           0.017 / 0.018
PXD004242                   0.936          0.926              0.936           0.010 / 0.001
cross-cohort: tune on PXD001064 -> 0.20, gives J 0.950 on PXD004242 (its optimum 0.957 at 0.225)
              tune on PXD004242 -> 0.225, gives J 0.601 on PXD001064 (its optimum 0.609 at 0.20)
```

At B's tuned margin, sensitivity and false-alarm rate are indistinguishable
from those at the 0.20 default (PXD001064 0.606/0.125 vs 0.609/0.124). On
the repo's own planter, the J-optimal margin is 0.175 for PXD001064 and
0.225 for PXD004242. That is consistent with "0.20 is near the best point
on both", but not literally the best for either.

**How much it moves.** About 1–2 points of J. That is an order of magnitude
smaller than the planting artefact in F2.

**Fix.** One sentence in the paper: "The default margin was chosen on these
two cohorts; nested donor-split tuning puts the resulting optimism at about
0.01–0.02 in Youden J, and a margin tuned on either cohort is within 0.008 J
of the other's optimum."

### F4. Figures 2D and 4B show an assignment better than every random draw

**What is wrong.** The "timepoint-rotation" assignment drawn in Figure 2D,
Figure 4B and the pair-score export is described as deterministic, which it
is. But it is not representative.

```text
draws with corrected AUC >= rotation 0.9839: 0 of 200 | max 0.983
draws with sens99 >= rotation 0.7760: 0
```

(from `audit/outputs`; one-liner over `pxd004242_draws.tsv` written by `a2`)

**How much it moves.** Figure 2D prints 0.9839 in bold and the
representative 0.9733 in small grey. A reader takes away the best case. The
caption's "AUC rises … from 0.8425 to 0.9839" is that best case.

**Fix.** Draw Figures 2C/D and 4B from the random assignment whose
corrected AUC is the median of the 200, and caption it as such. Or keep the
rotation and say in the caption that it lies above all 200 random
assignments. Then `export_pxd004242_pairs.py` should freeze that draw's
AUCs instead of the rotation's.

### F5. The two headline AUCs use different pair definitions (audit question 3)

**What differs:**

| | PXD001064 0.9692 | PXD004242 0.9733 |
|---|---|---|
| Negatives | All cross-batch different-donor pairs, *including* 137 MZ/DZ twin pairs | Cross-plate pairs *within one plate group* only |
| Positives | 93 pairs, exactly one per donor | About 674 pairs from 43 donors (about 16 per donor), pair-weighted |
| Summary | One deterministic dedup | Median over 200 run assignments |

**Demonstration:**

- Without twins among the negatives, PXD001064 gives 0.9707,
  sensitivity 0.688 (`a1`).
- With all cross-plate negatives, PXD004242 gives median 0.9674
  [0.9573, 0.9753] and sensitivity 0.661 [0.615, 0.712] (`a2`).

The within-group restriction is not conservative. Per-plate centring makes
different donors who share a plate group anti-correlated (`a4`):

```text
different-donor cross-plate: within-group median=-0.0359  cross-group median=-0.0032
plate-pair-stratified AUC=0.9836  vs pooled 0.9839        (no plate-pair leak)
same-date negatives score lower (-0.068 vs -0.035); positives share a date 3.2% vs 1.5%  (no date leak)
```

**Audit question 5.** Restricting to within-group pairs does remove the
confound it was meant to remove. There is no leak at the plate-pair or
date level. But plate-group structure *does* survive into the score,
because within-group pairs sit about 0.03 lower. That structure is what
makes the restricted negatives easier: +0.006 AUC and +3.6 points of
sensitivity at 99% specificity compared with all cross-plate pairs. The
tool itself scores rivals against *all* other batches, so its operating
point is the unrestricted one.

**Fix.** Report a matched pair in the main text. Both use unrelated-only
negatives across all cross-batch / cross-plate pairs: **0.9707 vs 0.9674**,
with donor-cluster CIs. Keep the within-group number as a design-aware
secondary analysis. Replace "same answer" with "consistent within
uncertainty". The CIs (F8) are far wider than the difference.

### F6. Label permutation vs physical swap (audit question 4)

**When the permutation is faithful.** The scores depend only on batch,
never on donor labels, so a within-batch label permutation produces exactly
the data a reciprocal tube swap *after plating* would. For that error, the
simulation is neither optimistic nor pessimistic. (Swaps confined to the
first batch, F2, are a separate problem.)

**Error modes it does not cover** (`a6_robust.py`, shipped tool, margin
0.20):

```text
intruder (a run from a donor not in the cohort, labelled as an enrolled donor), 40 trials:
   PXD001064 {'review': 27, 'swap_suspected': 13}
   PXD004242 {'swap_suspected': 33, 'review': 6, 'ok': 1}
off-by-one label shift across one batch:
   PXD004242 batch 3C: 12/12 relabelled runs -> swap_suspected
   PXD001064 batch 130109: of 36 relabelled, 17 swap_suspected, 12 review, 4 no_reference, 3 ok
all labels scrambled: PXD001064 128 swap / 74 review / 2 ok;  PXD004242 253 swap / 47 review
```

**Intruders: optimistic.** In PXD004242, 83% of intruders are called
`swap_suspected` *and the tool names an enrolled donor as the source*, who
is necessarily wrong. A rival run is compared against about 290 others at a
per-pair 99% threshold, so some stranger clears it most of the time. The
verdict text and README ("another donor matches it better") present that
best-of-290 as an identification.

**Reciprocal swaps: slightly pessimistic.** They create collateral that a
one-sided mislabel would not.

**Fix.** State that the simulation models reciprocal within-plate swaps
only. Report the intruder result. Change the `swap_suspected` description so
it does not imply the named donor is the true source.

### F7. Preprocessing without held-out donors (audit question 2)

**Demonstration** (`a4_leak.py`). There are 50 donor splits; PXD001064 is
split by twin pair and PXD004242 is stratified by plate group. The
99%-coverage set, batch medians and z-scaling are fitted on training donors
only and applied to held-out donors. AUC is computed on held-out pairs:

```text
                     transductive  held-out fit  drop (median)
PXD001064              0.9675        0.9624        0.006
PXD004242 rotation     0.9839        0.9687        0.015
```

**Threshold transfer** (`a8_ci.py`). A 99% threshold from one half's
negatives, applied to the other half under the tool's own transductive
scoring, gives realised specificity **0.9896** in both cohorts. Sensitivity
is 0.612 (PXD001064) and 0.684 (PXD004242).

**Design constraint.** In PXD004242, if held-out runs are corrected with
plate medians from *other* donors, specificity collapses to 0.91. The plate
composition is donor-nested, so the medians no longer centre the held-out
runs. The method therefore requires batch medians estimated on the full
batch, including the runs under test. The tool already does this. The paper
should say it is a requirement, not an implementation detail.

**Verdict.** The unsupervised steps do not leak donor labels in any way
that matters. The decision threshold uses labels, but it transfers.

### F8. No confidence intervals anywhere

Every interval in the paper is an assignment or panel range, and the
captions say so, but these ranges do not estimate sampling uncertainty.
Donor-cluster bootstrap (`a8`, 300 resamples):

```text
PXD001064: AUC 0.9692 [0.9374, 0.9855]  sens99 0.677 [0.422, 0.803]
PXD004242: AUC 0.9731 [0.9488, 0.9871]  sens99 0.668 [0.515, 0.801]   (one random draw)
```

The operating point is uncertain by about ±0.15–0.20. Any comparison of
0.677 with 0.697, or a reading of Figure 4C as two different operating
points, is within noise. **Fix:** add cluster-bootstrap CIs to Table 1.

### F9. Stated limits that are understated

- **"Single-run donors cannot be assessed" understates the limit.** In the
  shipped PXD001064 run, 44 of the 45 `no_reference` runs belong to donors
  with two runs *on the same acquisition day*. The tool excludes
  within-batch comparisons, so those donors cannot be assessed either
  (`a6`: "donor has 1 run: 1, donor has >=2 runs all in ONE batch: 44").
  The tutorial says `no_reference` "is large here because many donors were
  sampled only once", which is false: 115 of 116 donors have two runs. The
  real limit is that **donors whose runs all share one batch cannot be
  assessed**. That is 19% of PXD001064's runs.
- **"Identical twins" understates the relatedness problem.** DZ twin pairs
  clear the tool's own 99% threshold 19.3% of the time, against 0.8% for
  unrelated pairs, roughly 24× (`a6`/inline check; MZ 46.3%, same donor
  67.7%). Any first-degree relative in a cohort (siblings, parent–child) is
  a likely source of false rivals. The tutorial mentions relatives; the
  README and the paper's limits list do not.
- **Rival adjustment off by default: confirmed as the right call.** At
  matched false-alarm rates the margin dominates (`a5`):
  - PXD004242: `--rival-adjust 10` gives sensitivity 0.61 at 0.5%
    uninvolved flags; margin 0.30 alone gives 0.91 at 0.3%.
  - PXD001064: adjust 10 gives 0.27 at 2.1%; margin 0.30 gives 0.66 at
    1.1%.

  "Halved sensitivity" understates the cost; in PXD001064 it is about
  two-thirds. Also, above about 10 the setting is silently a no-op: the
  "keep 20 pairs in the tail" cap binds, and adjust 100 and 1000 give
  identical output. The CLI help should say so.

### F10. Tool defects found by reading and running the code

All of these are reproduced in `a6_robust.py`:

| Defect | Demonstration | Fix |
|---|---|---|
| **A failed run is reported as `no_reference`.** If 95% of a run's values are missing, every pair falls under `min_overlap`; `own` is empty, so the verdict reads "donor has no other run". | `failed run ...: verdict=no_reference n_ref=0` | Give it its own verdict (e.g. `unscorable`), driven by the run's finite-pair count. |
| **One bad run changes everyone's features.** The 99%-coverage selection is cohort-wide. | That same run dropped the panel from 273 to 259 features for all runs. | Report it in `summary.json`; consider excluding unscorable runs before selection. |
| **Tiny batches give silent artefacts.** A run alone in its batch is centred to zero: its own-donor score goes 0.467 → −0.027 and it becomes `review`. Two runs alone in a batch become anti-correlated (r = −0.941), and one of them is called `swap_suspected`. | `two runs alone in a batch: r between them=-0.941; verdicts=['swap_suspected','ok']` | Warn or refuse below a minimum batch size (at least 5–10 runs, say). |
| **The single-batch warning points the wrong way.** The CLI says scores "will be optimistic". Dropping batch correction actually multiplies false swap calls. | `no batch annotation: PXD001064 swap_suspected 67 (vs 7); PXD004242 22 (vs 5)` | Reword: expect many more false flags. |
| **Tutorial claims contradicted by the tool's own behaviour.** It says `review` "usually means the run itself is weak", but review runs barely differ from ok runs in missingness (0.205 vs 0.187) and not at all in total intensity (8.831 vs 8.839). It says a stranger match "is not something noise usually produces", but uninvolved clean runs are flagged at 1.5–4% and intruders get a stranger 33–83% of the time. | `a6`, F2, F6 | Remove both sentences, or support them with a measurement. Run quality was assessed only through these two proxies. |
| **Dedup depends on row order.** Duplicate donor/visit runs are resolved by the first run in matrix row order. | `keep="last"` gives AUC 0.9690, sensitivity 0.667 (vs 0.9692 / 0.677) | Document the rule; it does not change any conclusion. |

**Checked and fine:**

- Verdicts are invariant to row order ("row-order invariant: True").
- `pair_table` counts each pair once (`triu_indices`).
- NaN scores are excluded consistently.
- The CLI reproduces the tutorial output exactly (7 / 53 / 45 / 126,
  threshold 0.3161).
- With all labels scrambled, the tool flags almost everything, as it should.

### F11. Figure 6B overstates swap detection compared with the shipped tool

`analysis_pxd004242_swap_reference_curve.py::prepare_query` restricts the
candidate donors to the query's own plate group, a 1-of-about-12 closed set
that uses the plate layout. The tool searches all donors. "Label
disagreement" also has no clean-run counterpart. My re-implementation
(`a7_fig6b.py`, 5 draws × 200 within-plate swaps, all visits):

```text
                        in-group candidates   all 43 donors
swap recovery                 0.961               0.928
clean-run disagreement        0.037               0.073     <- the false-alarm rate Fig 6B omits
```

(Excluding the swap partner from references made no difference.)

Figure 6B's numbers are not the tool's: 0.974 at 7 visits against the
tool's 93–98%. **Fix:** caption it as closed-set identification within a
plate group, and add the clean-disagreement rate. Or redraw it with the
shipped tool and a visit-count subsample.

### F12. Smaller items

- **Figure 2 standardisation.** Standardising each panel by its own null is
  a per-panel affine map. It leaves the AUC unchanged, and the "gap in null
  SD" is a d′-type separation, so the standardisation is defensible and not
  in itself misleading. What misleads is the panel title (F1).
- **Figure 5B mixes estimands.** Random-273 is evaluated under the rotation
  assignment (the best case, F4); adaptive and top-273 use 200 random
  assignments. The bias works *against* the claim, so the conclusion holds.
  The correct like-for-like value is 0.9347/0.466 (`a9`).
- **Contaminant groups in PXD001064 but not PXD004242.** 51 of the 273
  PXD001064 features carry `CONT_` identifiers, 6 of them pure
  contaminants (e.g. porcine trypsin P00761, BSA P02769). PXD004242's
  extraction drops contaminants. Removing them changes almost nothing
  (`a4`: 0.9676 or 0.9685, sensitivity 0.677 unchanged). Disclose the
  inconsistency or harmonise it.
- **Cohort description.** "43 donors × 7 visits": participant 44 has 6.
- **Cosmetic.** In Figure 5B the "0.715" label collides with its error bar;
  in Figure 4A/B the subtitle text overlaps the histograms.

---

## 3. What the paper must say that it currently does not

1. **The 2×2 ablation.** Z-scaling and batch correction together give the
   gain. Z-scaling alone reaches 0.93 / 0.96 AUC. Batch correction carries
   the 99%-specificity end (0.41→0.68 and 0.64→0.78).
2. **Matched-estimand AUCs with donor-cluster CIs.** 0.9707
   [about 0.94–0.99] vs 0.9674. Say that the within-group PXD004242
   estimand is slightly favourable, and why.
3. **Uncertainty on the operating point.** Sensitivity at 99% specificity
   0.677 [0.42, 0.80] and about 0.67–0.70 [0.52, 0.80]. Say these are
   donor-cluster CIs, unlike the assignment ranges.
4. **Planting-design-independent swap validation.** Swaps spread across
   batches, with caught, uninvolved-flag and collateral-flag rates reported
   separately, and the clean-cohort `review` rate (17–23%).
5. **Margin-tuning optimism.** About 0.01–0.02 in Youden J (nested); the
   margin transfers across the two cohorts.
6. **What the swap model covers.** Reciprocal within-plate swaps only.
   Intruders are often attributed to an innocent enrolled donor.
7. **Limits stated at their real size:**
   - donors whose runs share one batch are unassessable (19% of PXD001064
     runs);
   - relatives beyond MZ twins (DZ 19% above threshold);
   - small batches (fewer than about 5 runs) are unsafe;
   - batch medians must be computed on the cohort being checked.
8. **Figures.** Figures 2D and 4B show an assignment above all 200 random
   draws; switch to a median draw or say so. Figure 6B is closed-set
   in-group identification, not the tool.

---

## 4. Analyses outside the audit scope

- **Not re-derived:**
  - the leave-one-acquisition-day-out (0.9681 / 0.9705) and
    leave-one-plate(-round)-out analyses;
  - the fixed-replicate diagnostics;
  - the 99.75% within-plate swap result (`analysis_pxd004242_swap_simulation.py`);
  - the swap reference curve's exact published values (only the comparison
    in F11 was reimplemented).

  None of these is among claims 1–7. Their scripts share preprocessing with
  the reproduced analyses, but their numerical results were not verified.
- **Sex check (PXD062541) and timestamp probe:** not looked at. They are
  presented as negative or supporting results, not claims.
- **Twin-pair labels vs the publication's Table S1:** not checked, because
  Table S1 is not in the repo. Only internal consistency was verified
  (zygosity agrees within all 58 pairs; 36 MZ, 22 DZ).
- **PXD004242 structure:** the SDRF confirmed 43 longitudinal donors, 300
  samples, 1,200 runs, plate groups of 12/13/12/6 donors and no donor in
  two groups. The tracked matrix was verified to have the same values as a
  fresh extraction from the zip. The 1,294 vs 1,290 file counts were not
  reconciled.
- **I/O readers:** the DIA-NN reader and the MaxQuant zip reader were not
  exercised on real DIA-NN or MaxQuant outputs. Only the generic matrix path
  was run end to end.
- **Not read:** `archive/`, the CI workflow, `plot_fig1_design.py` and
  Figure 1, `proteotrace_style.py`, `PUSH.md`.
- **Real-world swap rate:** no cohort with known real swaps was available.
  The intruder, off-by-one and scrambled-label cohorts in F6 are
  constructed, not observed.

---

## 5. Fixes applied (2026-10-02, branch `audit`)

| Finding | What changed | Where |
|---|---|---|
| F1 attribution | 2×2 ablation (raw / scaling / batch / both) in every headline; Figure 2 has three columns; README and log reworded | `analysis_primary_estimand.py`, `plot_fig2_batch_correction.py`, README |
| F2 tool validation | Swaps drawn over all batches; 30 repeats; caught / uninvolved / collateral reported separately; clean-cohort verdicts reported | `validate_tool.py`, `data/tool_validation.json`, README, tutorial §8 |
| F3 margin optimism | Number stated in README, CLI help and log | README, `proteotrace/cli.py` |
| F4 favourable assignment | Figures 2 and 4 use the median of 200 random draws; exporter checks against independently verified values | `export_pxd004242_pairs.py`, `data/PXD004242/PXD004242_primary_pairs_median.tsv` |
| F5 estimand | Primary estimand is the same in both cohorts (unrelated negatives, all cross-batch pairs); old numbers kept as secondary | `analysis_primary_estimand.py`, all figures |
| F6 swap model | README and tutorial state what the simulation covers; `swap_suspected` text no longer implies identification | README, tutorial §5, `proteotrace/report.py` |
| F7 held-out | Stated; batch medians must come from the cohort itself (core docstring, README, tutorial §7) | README, `proteotrace/core.py` |
| F8 CIs | Donor-cluster bootstrap CIs in README, log and Figure 4C | `data/*primary*.json`, `plot_fig4_ceiling.py` |
| F9 limits | Relatives, same-batch donors, small batches, outside samples in README / tutorial / `summary.json` caveats | README, tutorial §10, `proteotrace/report.py` |
| F10 tool defects | `unscorable` verdict; small-batch warning; single-batch warning corrected; `--rival-adjust` cap documented; tutorial claims corrected. Also fixed: the tutorial called two DZ twin rivals "identical twins" | `proteotrace/*`, tutorial §4–5 |
| F11 Figure 6B | Replaced by the shipped tool's caught / uninvolved / collateral rates by visits per donor | `analysis_tool_visit_curve.py`, `plot_fig6_qc_use.py` |
| F12 Figure 5B | All three panels evaluated on the same 200 draws and the same estimand | `plot_fig5_transfer.py` |

**Not changed:**

- Contaminant groups remain in the PXD001064 panel; removing them changes
  nothing measurable, so this is a disclosure item for the methods text.
- The dedup rule (first run in matrix order) is unchanged and documented.
- The secondary analyses that no main figure reads (leave-one-out folds,
  fixed-replicate diagnostics, sex check) were not re-run.

### Follow-up, later on 2026-10-02

These items came out of checking which external data exist. Details are in the
analysis log.

| Item | What changed |
|---|---|
| Third cohort | PXD062541 (another lab, Exploris 480 / DIA-NN), run through the tool at its defaults: AUC 0.892 [0.792, 0.956], sensitivity at 99% specificity 0.25. The method transfers, but clearly less well without batch information and at lower depth (`analysis_pxd062541.py`). |
| Decimal-comma bug | The DIA-NN reader turned `7,04E+12` into NaN and took the maximum of text values; fixed, with byte-identical results on the other two cohorts (0.3.0). |
| Coverage rule | A protein missing from a single run is always kept. Identical results from 100 runs upward; Figure 6B at two visits moves from 58.9% to 59.8%. |
| Prior work | Geyer et al. 2021 (MCP) already re-identified the PXD004242 participants from plasma proteomes and named sample mix-up as a use. This is cited in the README; claims of novelty should rest on the batch handling, the calibrated threshold with measured false alarms, the cross-cohort test and the tool. |
