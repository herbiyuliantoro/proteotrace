# Extension analysis plan: individuality weighting and genotype-like peptides

Written and committed **before** any of the analyses below were run, so that
the success criteria cannot be adjusted after seeing results. Deviations from
this plan will be reported as deviations.

## Rules that apply to both ideas

1. **Development cohorts:** PXD001064 and PXD004242. Every choice (weights,
   thresholds, peptide filters, how scores are combined) is fixed using these
   two cohorts only.
2. **Held-out test cohort:** PXD062541. It is evaluated **once**, with the
   choices frozen after development, and the result is reported whatever it is.
3. **Labels never enter the score for the pairs being evaluated.** Any quantity
   estimated from donor labels (variance components, peptide selection) is
   estimated on training donors or on a different cohort, never on the donors
   whose pairs are scored.
4. **Primary metric:** sensitivity at 99% pairwise specificity on the primary
   estimand (same-donor vs unrelated cross-batch pairs; all pairs in
   PXD062541). Secondary: AUC. Comparisons with the current score are made on
   identical pairs, with donor-cluster paired bootstrap 95% CIs (1,000
   resamples) for the difference.
5. **Success criterion (per idea):** in development, the improvement in the
   primary metric has a paired 95% CI excluding zero in at least one
   development cohort and a point estimate ≥ 0 in the other. If met, the
   method is evaluated once on PXD062541.
6. **The released tool's defaults are not changed** for the current
   manuscript. Results are reported as an extension.

## Idea 2: weighting proteins by individuality

- Per protein, on the batch-corrected, z-scaled values used by the current
  score, estimate the intraclass correlation ICC = σ²between / (σ²between +
  σ²within) from a one-way random-effects decomposition over donors with ≥ 2
  runs. Negative estimates are set to 0.
- **Weighted score:** weighted Pearson correlation of two runs, with protein
  weights equal to the ICC, over proteins measured in both runs.
- **Design A, within-cohort cross-fitting (development):** split donors into two
  halves (twin pairs kept together in PXD001064; stratified by plate group in
  PXD004242), estimate ICC on one half and score pairs of the other half, then
  swap. Repeat for 20 random splits; report the median difference and the
  paired bootstrap CI from the first split. PXD004242 uses the assignment
  closest to the median (seed 20261090).
- **Design B, transfer without labels in the target (the deployable version):**
  ICC is estimated on PXD004242 (all donors; 6–7 visits each) and applied to
  PXD001064 (development) and PXD062541 (test) by matching UniProt
  accessions (first accession of each protein group). Proteins without a
  transferred weight receive the median transferred weight.

## Idea 3: genotype-like peptides

- **Feature:** for each peptide, log2(peptide intensity) − log2(intensity of
  its protein group) per run, i.e. the peptide-to-protein ratio. Peptides
  must be observed in ≥ 80% of runs of the cohort.
- **Selection by individuality:** ICC of the ratio, estimated as in idea 2 on
  training donors (design A) or on PXD004242 (design B, matched by peptide
  sequence).
- **Peptide score:** ICC-weighted Pearson correlation over the ratio features
  of two runs, after the same batch-median subtraction and z-scaling as the
  protein score.
- **Combined score:** mean of the protein score and the peptide score.
- **Biological check (not used for selection):** the fraction of the 50
  highest-ICC peptides that overlap a UniProt natural-variant position,
  compared with the fraction among all eligible peptides.
- **Data:** PXD004242 `peptides.txt` (MaxQuant), PXD062541 DIA-NN precursor
  report, PXD001064 `OpenSWATHResult_feature_alignment_requant.tsv` (PRIDE).
- **Out of scope:** re-searching raw files against a variant protein database
  (hundreds of GB of raw data); noted as future work.

## Expected limits (stated in advance)

- Neither idea can separate monozygotic twins, who share their genome; idea 3
  should, if anything, make MZ pairs *more* similar to same-donor pairs.
- Gains in two-visit cohorts are bounded by having one reference pair per run.

---

## Results (added after the analyses; the plan above is unchanged)

### Deviations from the plan

- **Protein reference for peptide ratios:** the median log2 intensity of the
  protein's peptides in each run, rather than the protein-group intensity of
  each search engine, so that MaxQuant, DIA-NN and OpenSWATH use one definition.
- **OpenSWATH filtering (PXD001064):** decoys removed; features with
  m_score ≤ 0.01 kept, together with requantified features (which carry no
  m_score).
- **Exploratory, not pre-registered:** a combined score using the
  idea-2-weighted protein score, and the relatedness ladder for peptide scores.

### Idea 2: individuality weighting (`analysis_individuality.py`)

Success criterion met. Sensitivity at 99% specificity, weighted minus current
score:

| Cohort | Design | Difference | Paired 95% CI |
|---|---|---|---|
| PXD004242 | A, split 1 (median of 20 splits) | +0.140 (+0.133) | +0.068 to +0.212 |
| PXD001064 | A, split 1 (median of 20 splits) | +0.097 (+0.097) | −0.006 to +0.264 |
| PXD001064 | B, weights from PXD004242 | +0.043 | −0.024 to +0.135 |
| **PXD062541 (held-out test)** | B, 59/133 proteins matched | **+0.071** | −0.179 to +0.224 |

All 20 splits were positive in both development cohorts. The held-out test
points in the same direction but is not significant; only 59 of 133 proteins
had a transferred weight.

### Idea 3: genotype-like peptides (`analysis_genotype_peptides.py`)

Success criterion **not met**, so the method was not evaluated on PXD062541.
Combined score minus current score, median of 20 splits: −0.071 (PXD004242)
and −0.065 (PXD001064). The peptide score alone was much weaker (−0.537 and
−0.414).

The following observations are consistent with a genetic contribution to
peptide similarity:

- Of the 50 peptides with the highest ICC in PXD004242, 24 (48%) overlap a
  UniProt natural-variant position, against 255 of 1,176 (22%) other eligible
  peptides (one-sided Fisher p = 5.4 × 10⁻⁵).
- In PXD001064 (design B weights), the median peptide score was +0.140 for MZ
  twins, +0.121 for the same donor, +0.058 for DZ twins and −0.009 for
  unrelated pairs. MZ twins and repeat samples have similar median peptide
  scores, while the DZ-twin median is lower. The protein-score medians are
  +0.381 for the same donor and +0.293 for MZ twins. These patterns and the
  variant-position overlap do not by themselves establish a genetic cause;
  donor genotypes were not measured in this reanalysis.
