# 🖼️ Figures

What each of the six main figures and the supplementary Figure S1 claims, how
it is built, and the draft captions. The design rules are here so a replacement figure still matches the
rest of the set.

---

## Main figures

Six figures, one claim each. Nothing in a main figure is set below 7.5 pt: the
figures carry marks, axis labels, panel titles and a few direct values, and
every design fact, interval definition and sample size lives in the caption
below. The earlier four-figure set is kept under `archive/`.

```sh
cd proteotrace
for f in plot_fig*.py; do python "$f"; done
```

### Draft captions

Conventions used on the figures. Cohorts are named by PRIDE accession:
PXD001064 (twin cohort), PXD004242 (weight-loss cohort) and PXD062541
(epicatechin cohort, Figure S1). Positive pairs are two runs from the same
donor; negative pairs are runs from different donors. *Specificity* is the
share of negative pairs below a threshold, *sensitivity* the share of
positive pairs above it, and *AUC* the probability that a random positive
pair scores higher than a random negative pair. MZ, monozygotic (identical)
twins; DZ, dizygotic (fraternal) twins. Panel titles state the finding of
each panel.

**Figure 1. How the similarity score is computed, and which pairs it is read on.**
(A) The three steps (protein quantification, batch-median correction,
pairwise correlation), shown on a real slice of PXD001064:
(1) fourteen proteins measured across ten runs from two acquisition days;
(2) the same values after the per-batch median of each protein is subtracted;
(3) the two runs of one donor, each protein scaled to unit variance across the
cohort, correlated over all 273 proteins.
Rows are ordered by their between-batch difference, so the block structure the
correction removes is legible. The colour scale is clipped at the 80th
percentile of absolute values. The example pair is the same-donor cross-batch
pair whose score is closest to the median of that class, r = 0.381. Donor
labels are used only to evaluate the score, never to compute it.
(B) The pairing rule. Only pairs on different batches (acquisition days in
PXD001064, plates in PXD004242) are read, across plate groups as well as
within them. This is because the tool compares each run with every other
batch. Same-batch pairs are excluded. In PXD001064 the 137 twin pairs are left
out of the primary negatives and analysed in Figure 3.

**Figure 2. What each preprocessing step adds.**
Pair scores in PXD001064 (top; same-donor vs unrelated cross-batch pairs,
273 protein groups) and PXD004242 (bottom; cross-plate pairs, 150 protein
groups). The PXD004242 assignment is the median of 200 random one-run-per-sample
assignments.
(A) Raw: log2, centred on each run's median. (B) Each protein additionally
scaled to unit variance. (C) Additionally, the per-batch median of each protein
subtracted.
Each panel is standardised against its own different-donor distribution,
because raw correlations all lie near 0.95 while corrected scores span −0.5 to
0.7. This per-panel affine map leaves the AUCs unchanged.
AUC rises 0.7836 → 0.9329 → 0.9707 (PXD001064) and 0.8389 → 0.9588 → 0.9674
(PXD004242). Scaling carries most of the AUC, but batch correction carries
the strict end: sensitivity at 99% specificity rises from 0.41 to 0.69, and
from 0.59 to 0.66 (median of 200 assignments; Table 1).

**Figure 3. Similarity declines with genetic relatedness.**
Corrected cross-batch scores in PXD001064 for four relations: same donor
(n = 93), monozygotic twin (n = 80), dizygotic twin (n = 57) and unrelated
(n = 23,253) pairs. Points are medians, thick bars the 25th–75th and thin bars
the 5th–95th percentiles; printed values are the medians. The same-donor and
MZ distributions overlap, so the
screen does not separate identical twins. DZ twins also clear the 99%
threshold in 19% of pairs, against 0.8% for unrelated pairs. These are stated
limits of the method rather than a tuning problem.

**Figure 4. A strict threshold costs sensitivity.**
(A, B) Corrected score distributions on the primary estimand, in PXD001064 (A)
and the median PXD004242 assignment (B). Thresholds sit at the 95th and 99th
percentile of different-donor pairs; the same-donor mass left of each line is
what a screen at that specificity misses.
(C) The resulting operating points. PXD004242 points are medians over 200
assignments. Bars are donor-cluster bootstrap 95% CIs (1,000 resamples; twin
pair is the cluster in PXD001064, donor in PXD004242). The PXD004242 CI is
computed on the median assignment.

**Figure 5. The external result does not rest on arbitrary choices.**
(A) 200 random assignments of one run per biological sample in PXD004242,
through each preprocessing step. Grey lines join the same assignment. The
dashed line is the PXD001064 estimate of 0.9707.
(B) Three protein panels, all evaluated on the same 200 assignments:
- adaptive panels (proteins observed in ≥99% of runs; median 147);
- the 273 highest-coverage proteins;
- 273 randomly chosen proteins (twenty panels, cycled across assignments;
  median coverage 25%).

High-coverage proteins carry the signal; arbitrary, mostly sparse panels
lose it. Bars are 2.5–97.5% assignment (and, for random-273, panel) ranges.
They are not confidence intervals.

**Figure 6. What the screen delivers in use.**
(A) Donor discrimination against the gap between visits, with the pre-to-post
and post-only contrasts shown separately. Medians over 200 assignments; the
dashed line is the PXD001064 estimate of 0.9707. Identity is recovered across
the 12% weight-loss intervention and out to 60 weeks. Each gap rests on about
32 same-donor pairs, so the per-gap differences are within noise; points are
not joined, because the gaps are irregularly spaced.
(B) The shipped tool (`proteotrace check`, margin 0.20) on planted reciprocal
same-plate swaps. Each donor is subsampled to the given number of visits, with
12 swaps per repeat and 20 repeats per point. Three rates are shown:
- planted runs called `swap_suspected` ("swapped runs", i.e. sensitivity);
- runs of donors not involved in any swap ("unaffected runs", i.e. the
  false-positive rate);
- the swapped donors' other runs ("other runs of swapped donors"), which are
  compared against a wrong reference by construction.

Bars are 2.5–97.5% ranges over repeats. These are simulated swaps, not
validation on a cohort of unknown status.

### Supplementary figure

**Figure S1. The third cohort (PXD062541), run through the tool at its defaults.**
PXD062541 comes from a third laboratory (Orbitrap Exploris 480, DIA-NN): 28
adults, plasma before and after epicatechin supplementation. No batch
information is deposited, so there is no batch-correction step. It is also
shallower than the other two cohorts (median 44% of values missing per run).
(A, B) Same-person (n = 28) and different-people (n = 1,512) pair scores
without correction (A) and with each protein scaled (B, the tool's default),
standardised as in Figure 2.
(C) Sensitivity analysis of the protein-selection rule. AUC and same-person
pairs found at the 99% cut-off when proteins must be seen in at least 99%
(* tool default), 95%, 90% or 80% of runs; the number under each threshold is
how many proteins were kept. Bars on the default are donor-cluster bootstrap
95% CIs. The default is the pre-specified primary result; the looser rules are
reported only as a sensitivity analysis (`analysis_pxd062541.py`).

### Two deliberate choices about what is plotted

- **Figure 2 standardises each panel against its own null.** On one raw axis the
  before-correction panel renders empty, because raw correlations occupy a
  0.87-0.98 band while corrected scores span -0.5 to 0.7. In null-SD units the
  reader sees the quantity each step changes. The map is affine per panel, so
  every printed AUC is the AUC of the untransformed scores.
- **Figures 2 and 4 draw the median assignment, not a favourable one.** The
  earlier timepoint-rotation assignment scored above all 200 random draws
  (AUC 0.9839 against a maximum of 0.983), so drawing from it showed a best
  case. It is kept in the analysis log as a sensitivity analysis only.
- **Figure 4 replaces the ROC curve.** A curve hugging the corner says only
  "good". Drawing the thresholds on the two score distributions shows *why*
  sensitivity stops near two in three, which is the number a user has to plan around.

Moved to supplementary: the ROC curve, the acquisition-day and plate-count
summaries, leave-one-plate-round-out, the sex check, the technical-replicate
diagnostic, and the secondary estimands (PXD001064 with twins among the
negatives, 0.9692; PXD004242 restricted to one plate group, 0.9733; the
timepoint-rotation assignment, 0.9839). Cohort descriptions and the full AUC and sensitivity numbers belong
in Table 1, where they can be read precisely.

Specification, so a replacement figure stays consistent:

- **Size.** Journal of Proteome Research (ACS) artwork rules: double-column
  figures 176 mm (499 pt; JPR allows 300–504 pt), Figure 3 single column 84 mm
  (238 pt; JPR allows up to 240 pt), all well under the 660 pt depth limit.
  `proteotrace_style.save()` refuses to write a figure that is too wide, has
  text under 4.5 pt or a drawn line under 0.5 pt, so the rules cannot drift.
- **Type.** Arial at 8 pt base, 9 pt bold panel titles, 11 pt bold panel
  letters. Nothing is set below 7.5 pt. Arial is registered from the Windows
  font directory when running under WSL and falls back to DejaVu Sans
  elsewhere, so the figures build on any platform; only the typeface changes.
  Headers are placed with `S.head()`, which offsets in typographic points
  rather than axes fractions, so the header block keeps its spacing at any
  panel height.
- **Colour.** Slots 1 to 3 of the validated reference palette: blue `#2a78d6`,
  orange `#eb6834`, aqua `#1baf7a`. That three-hue subset clears the all-pairs
  colour-vision gates (worst CVD dE 9.2, normal-vision dE 24.0). Greys are never
  a fourth category; they carry "raw", "baseline" or "unrelated". The plate
  heat map uses one hue light to dark, not a rainbow.
- **Marks.** No grid. Black left and bottom axes with outward ticks; density
  panels drop the left axis, since their height has no scale worth reading.
  Text is in ink, not in series colours; reference lines are dark grey and
  labelled. Direct value labels rather than a number on every point.
- **Titles and terms.** Each panel carries a short title stating its finding
  in plain scientific English; axes and legends use standard terms
  (sensitivity, specificity, AUC, Pearson r, MZ/DZ). Cohorts are named by
  PRIDE accession. Sample sizes, interval definitions and design facts live
  in the caption.
- **Export.** Vector PDF with TrueType (type 42) fonts for editable text,
  600 dpi TIFF (LZW) for submission (JPR asks for at least 300 dpi colour
  art), and 600 dpi PNG for viewing. TIFFs are not tracked in git; every
  `plot_fig*.py` rebuilds them.

Intervals in Figures 5 and 6 are run-assignment, panel-selection or resampling
ranges, not confidence intervals. Figure 4C is the exception: its bars are
donor-cluster bootstrap 95% CIs. Each caption says which is which.

## PXD004242 pair-score export

Figures 2 and 4 need the score of every pair, not a summary AUC, so
`export_pxd004242_pairs.py` re-runs the pipeline for the median random
assignment of `analysis_primary_estimand.py` (seed 20261090). It writes each
primary pair (cross-plate, any plate group) with its raw, scaling-only and
corrected score:

```sh
cd proteotrace
python analysis_primary_estimand.py     # writes the 200-draw summary and picks the median draw
python export_pxd004242_pairs.py
```

```text
       raw AUC 0.8389  (expected 0.8389)
    scaled AUC 0.9588  (expected 0.9588)
 corrected AUC 0.9674  (expected 0.9674)
seed=20261090 features=150 runs=300 pairs=41912 same_donor=666
```

The expected values were verified with the audit's independent implementation
(`audit/audit_lib.py`, which reads only the MaxQuant zip and the SDRF). The
script exits with an error if any AUC drifts by more than 5e-4, or if any count
changes. Output: `data/PXD004242/PXD004242_primary_pairs_median.tsv`.
