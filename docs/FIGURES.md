# 🖼️ Figures

What each of the six main figures claims, how it is built, and the draft
captions. The design rules are here so a replacement figure still matches the
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

**Figure 1. How the similarity score is computed, and which pairs it is read on.**
(A) The three steps, shown on a real slice of PXD001064: fourteen proteins
measured across ten runs from two acquisition days (1), the same values after
the per-batch median of each protein is subtracted (2), and the two runs of one
donor correlated across all 273 proteins (3). Rows are ordered by their
between-batch difference so the block structure the correction removes is
legible; the colour scale is clipped at the 80th percentile of absolute values.
The example pair is the same-donor cross-batch pair whose score is closest to
the median of that class, r = 0.381. Donor labels are used only to evaluate the
score, never to compute it. (B) The pairing rule. In PXD004242 each donor sits
in one plate group, so the primary comparison uses pairs from different plates
within a single group; pairs on the same plate, or spanning plate groups, are
excluded.

**Figure 2. Batch correction is the method, not a cleaning step.**
Pair scores before (A, C) and after (B, D) subtracting the per-batch protein
median, in PXD001064 (A, B; cross-batch pairs, 273 protein groups) and
PXD004242 (C, D; primary within-group cross-plate pairs, 148 protein groups,
timepoint-rotation assignment). Each panel is standardised against its own
different-donor distribution, because raw correlations all lie near 0.95 while
corrected scores span -0.5 to 0.7. The median gap widens from 0.97 to 2.94 null
standard deviations in PXD001064 and from 1.20 to 3.18 in PXD004242; AUC rises
from 0.7827 to 0.9692 and from 0.8425 to 0.9839.

**Figure 3. Similarity declines with genetic relatedness.**
Corrected cross-batch scores for same-donor (n = 93), monozygotic twin (n = 80),
dizygotic twin (n = 57) and unrelated (n = 23,253) pairs in PXD001064. Points are
medians, thick bars the 25th-75th and thin bars the 5th-95th percentiles. The
same-donor and MZ distributions overlap: the screen does not separate identical
twins, and this is a stated limit of the method rather than a tuning problem.

**Figure 4. A strict threshold costs sensitivity.**
Corrected score distributions in PXD001064 (A) and PXD004242 (B), with
thresholds placed at the 95th and 99th percentile of different-donor pairs; the
same-donor mass left of each line is what a screen at that specificity misses.
(C) The resulting operating points. PXD004242 panels use the timepoint-rotation
assignment, while its operating points in C are medians over 200 assignments
with the 2.5-97.5% assignment range as bars, not confidence intervals.

**Figure 5. The external result does not rest on arbitrary choices.**
(A) 200 random assignments of one run per biological sample in PXD004242, raw
against batch-corrected; the dashed line is the PXD001064 estimate of 0.9692.
(B) Adaptive panels (median 147 proteins), the 273 highest-coverage proteins and
273 randomly chosen proteins. Coverage quality, not panel size, carries the
signal. Bars are 2.5-97.5% assignment or panel-selection ranges.

**Figure 6. What the screen delivers in use.**
(A) Donor discrimination against the gap between visits, with the pre-to-post
and post-only contrasts shown separately; the dashed line is the PXD001064
estimate of 0.9692. Identity is recovered across the 12% weight-loss
intervention and out to 60 weeks. (B) Simulated within-plate swaps, scored against the number of
reference visits held for the donor. Label disagreement is near-certain once
several visits exist; correct donor recovery rises from 0.882 at two visits to
0.974 at seven. These are simulated swaps with known reference donors, not
validation on a cohort of unknown status. Bars are 2.5-97.5% resampling ranges.

### Two deliberate choices about what is plotted

- **Figure 2 standardises each panel against its own null.** On one raw axis the
  before-correction panel renders empty, because raw correlations occupy a
  0.87-0.98 band while corrected scores span -0.5 to 0.7. In null-SD units the
  reader sees the quantity the method actually changes.
- **Figure 4 replaces the ROC curve.** A curve hugging the corner says only
  "good". Drawing the thresholds on the two score distributions shows *why*
  sensitivity stops near 0.68, which is the number a user has to plan around.

Moved to supplementary: the ROC curve, the acquisition-day and plate-count
summaries, leave-one-plate-round-out, the sex check, and the technical-replicate
diagnostic. Cohort descriptions and the full AUC and sensitivity numbers belong
in Table 1, where they can be read precisely.

Specification, so a replacement figure stays consistent:

- **Size.** 180 mm wide (double column); the PDFs carry a 510.24 pt media box.
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
- **Marks.** Thin rules, no chart borders beyond left and bottom, recessive grid
  at `#e6e5e1`, direct value labels rather than a number on every point.
- **No embedded figure titles or footnotes.** Captions belong in the manuscript;
  every panel instead carries a one-line subtitle giving n and the design fact
  needed to read it.
- **Export.** 600 dpi PNG plus vector PDF with TrueType (type 42) fonts embedded,
  which is what JPR and MCP require for editable text.

Every interval drawn in these figures is a run-assignment, panel-selection or
resampling range. None of them is a statistical confidence interval, and the
panel subtitles say so.

## PXD004242 pair-score export

Figures 2 and 4 need the score of every pair, not a summary AUC, so
`export_pxd004242_pairs.py` re-runs the prespecified pipeline for the
deterministic timepoint-modulo-4 assignment and writes each primary pair with
its raw and corrected score:

```sh
cd proteotrace
python export_pxd004242_pairs.py
```

```text
       raw AUC 0.8425  (frozen 0.8425)
 corrected AUC 0.9839  (frozen 0.9839)
features=148 runs=300 pairs=8826 same_donor=759
```

The script recomputes both AUCs and exits with an error if they drift from the
frozen rotation result by more than 5e-4, so a figure can never be drawn from a
run that no longer reproduces the recorded analysis. Output:
`data/PXD004242/PXD004242_primary_pairs_rotation.tsv`.
