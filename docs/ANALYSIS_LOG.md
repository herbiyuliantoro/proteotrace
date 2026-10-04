# 📓 Analysis log

The running record: what was run, what it produced, and which conclusions were
retired on the way. Kept as written so the reasoning behind the frozen numbers
stays auditable.

For the headline results see the [main README](../README.md); for the figures
see [FIGURES.md](FIGURES.md).

---

## Corrections after the pre-preprint audit (2026-10-02)

An adversarial audit ([AUDIT_BRIEF.md](AUDIT_BRIEF.md),
[AUDIT_REPORT.md](AUDIT_REPORT.md)) re-derived every headline number from the
raw inputs with independent code (`audit/`). The central discrimination result
reproduced exactly. Five things changed.

**1. Attribution.** The "raw" baseline left out both the per-batch median and
the per-protein z-scaling, so the raw-to-corrected gain had been credited to
batch correction alone. `analysis_primary_estimand.py` now reports the 2×2
ablation (AUC / sensitivity at 99% specificity):

```text
                 raw            scaling only     batch only       both
PXD001064   0.7836 / 0.118   0.9329 / 0.409   0.8828 / 0.226   0.9707 / 0.688
PXD004242   0.8438 / 0.281   0.9579 / 0.591   0.9444 / 0.522   0.9674 / 0.661   (median of 200)
```

The paper claim becomes: per-protein scaling and within-batch median centring,
both unsupervised, together make routine plasma proteomes usable for
same-donor detection. Scaling carries most of the AUC; batch centring carries
the 99%-specificity end.

**2. One estimand for both cohorts.** The primary estimand is now:

- positives: same-donor pairs on different batches;
- negatives: *unrelated* different-donor pairs on different batches.

In PXD001064 the 137 MZ/DZ twin pairs are excluded from the negatives; they
are analysed in the ladder. In PXD004242, pairs across plate groups are kept,
because the tool scores rivals against every other batch. The within-group
restriction turned out to be mildly *favourable*: per-plate centring makes
different donors of one plate group anti-correlated, with median score −0.036
against −0.003 across groups.

```text
primary   PXD001064 0.9707 [CI 0.942, 0.989]  sens99 0.688 [0.508, 0.842]   93 / 23,253 pairs
          PXD004242 0.9674 (assignments 0.9573-0.9753)  sens99 0.661 (0.615-0.712)
                    median draw seed 20261090: CI AUC [0.939, 0.984], sens99 [0.538, 0.770]
secondary PXD001064 incl. twin negatives 0.9692 / 0.677
          PXD004242 within plate group   0.9733 / 0.697 (as before)
```

CIs are donor-cluster bootstrap (1,000 resamples; twin pair or donor is the
cluster).

**3. Figures use the median assignment.** The timepoint-rotation assignment
(0.9839) scored above all 200 random draws, so figures drawn from it showed a
best case. Figures 2 and 4 now use the median draw. The export is checked
against values verified independently.

**4. The tool validation was rebuilt.** The old planter filled its 12 swaps
from the alphabetically first batch in every repeat. It also counted the
swapped donors' *other* runs, which are flagged by construction, as false
alarms. `validate_tool.py` now draws swaps over all batches (30 repeats) and
reports caught / uninvolved / collateral separately. At margin 0.20:

```text
PXD001064  caught 69.8%  uninvolved 3.8%  collateral 76%   (old table: 80.6% / 12.2%)
PXD004242  caught 97.6%  uninvolved 1.8%  collateral 4.7%  (old table: 93.3% / 4.9%)
```

Figure 6B now shows the shipped tool by visits per donor
(`analysis_tool_visit_curve.py`): 58.9% caught at 2 visits, 97.5% at 7. The
earlier reference curve (88.2% → 97.4%) scored a closed-set identification
within one plate group, which is not what the tool does.

**5. Margin selection and held-out donors, quantified.**

- Margin selection: tuning on half the donors and testing on the other half
  costs 0.017 (PXD001064) and 0.010 (PXD004242) in Youden J.
- Cross-cohort transfer: a margin tuned on either cohort is within 0.008 J of
  the other cohort's optimum (`audit/a5_margin.py`).
- Held-out donors: fitting coverage, batch medians and scaling on half the
  donors costs 0.006 / 0.015 AUC on the held-out half (`audit/a4_leak.py`).
- Threshold transfer: a 99% threshold from other donors' pairs gives realised
  specificity 0.990 in both cohorts (`audit/a8_ci.py`).

The tool also gained an `unscorable` verdict for failed runs, a warning for
batches under 5 runs, a corrected single-batch warning, and accurate
descriptions of `no_reference` and `review` (version 0.2.0).

Retired by this entry: "batch correction is the method" as worded below, the
0.9692-vs-0.9733 comparison as a like-for-like replication, the 93.3% / 80.6%
swap table, and Figure 6B's reference-curve numbers as a description of the
tool.

## Third cohort, and two reader/selection fixes (2026-10-02, later)

**Third cohort: PXD062541.** This is the epicatechin study from the Otton lab
(Cruzeiro do Sul University): Orbitrap Exploris 480, DIA-NN, 28 obese adults,
plasma before and after supplementation. It was already in the repository for
the sex check.

The reporting rules were decided before any result was looked at:

- the primary result is the tool at its defaults;
- other coverage thresholds are a sensitivity analysis only;
- single-batch mode, because no batch information is deposited;
- `WB_1` is read as W1_B (see below);
- 6 planted pairs per repeat.

`analysis_pxd062541.py` writes `data/PXD062541_primary.json` and
`examples/PXD062541_samples.tsv`.

```text
primary (tool defaults, 133 proteins, 28 same-donor / 1512 different-donor pairs)
  raw AUC 0.8042  ->  default AUC 0.8919 [CI 0.792, 0.956]
  sensitivity at 95% / 99% specificity: 0.464 / 0.250 [CI 0.07, 0.63]
  clean cohort: ok 14, review 38, swap_suspected 4 (W12_B, W2_A, W6_A, WB_1)
  planted swaps (30 x 6 pairs): caught 37.8%, uninvolved runs flagged 3.5%,
                                other runs of swapped donors flagged 39.4%
sensitivity analysis (same pipeline, other --min-coverage)
  0.95: 180 proteins, AUC 0.9194, sens99 0.357
  0.90: 278 proteins, AUC 0.9355, sens99 0.500
  0.80: 428 proteins, AUC 0.9360, sens99 0.357
```

The method transfers, but clearly less well. Three conditions differ from the
two validation cohorts at once, so their effects cannot be separated here:

- there is no batch information, so no correction;
- the data are shallower (median 44% of values missing per run);
- there are only two visits per donor.

**WB_1.** The DIA-NN report names one run `WB_1`. The SDRF has no `WB_1`, and
has a `W1_B.raw` that is otherwise absent from the report. These are the only
mismatches in either direction, so `WB_1` is read as W1_B. This is an
inference from file names, not a confirmation from the submitters.

WB_1 is among the four flagged runs. Its score to W1_A is 0.019, and its
nearest runs are W6_B (0.277) and W2_A (0.272). At this cohort's depth those
values are within the range of unrelated pairs, so the data cannot tell a weak
signal from a genuine mislabel. Only the submitters could settle it.

**Fix 1, decimal commas (proteotrace 0.3.0).** The report writes `7,04E+12`.
The DIA-NN reader turned every value into NaN, then stopped with a misleading
"only 0 features" error. Before the pivot, it also took the maximum of *text*
values.

The readers now parse decimal commas when a column never mixes "." and ",".
They convert before pivoting, and they stop with a clear error if no positive
number could be read.

Checks:

- the generic reader on PXD001064 and the MaxQuant reader on the PXD004242
  zip give byte-identical matrices to before;
- the DIA-NN matrix equals an independent manual conversion.

**Fix 2, coverage rule for small cohorts (0.3.0).** "Observed in ≥99% of runs"
demands a complete protein below 100 runs, so one poor run (M7_B, 85% missing)
cut the PXD062541 panel to 81 proteins (133 under the new rule). A protein missing from a single
run is now always kept. From 100 runs upward this is identical to the old rule:
`validate_tool.py` reproduces `data/tool_validation.json` byte for byte. The
Figure 6B visit curve changes only at two visits (86 runs), from 58.9% to 59.8%
caught; three to seven visits are unchanged.

**Fix 3, panel smaller than `--min-overlap` (0.3.1).** The coverage step only
refused panels below 20 proteins, while a pair needs `--min-overlap` (30) shared
proteins to be scored. A panel of 20–29 proteins therefore passed the first
check, left every pair unscored, and stopped with a misleading "no
different-donor pairs" error. `check` now stops right after the coverage step
with a message naming both settings. The 20-protein floor stays for users who
lower `--min-overlap`, so the effective minimum is max(20, `--min-overlap`).

Checks:

- a synthetic 25-protein cohort stops with the new message at the default and
  runs with `--min-overlap 20`;
- PXD001064 gives the same verdicts as before (126 consistent, 53 review,
  7 possible swaps, 45 without reference).

## Earlier framing (retired 2026-10-02, kept for the record)

The current paper should make one focused claim: an unsupervised, acquisition-day
normalization makes routine plasma proteomes useful for detecting whether two
runs are from the same donor. In the PXD001064 pilot, the cross-batch AUC was
**0.7827 before correction and 0.9692 after correction**. The first number must
remain in the paper: batch correction is the method's main contribution, not a
cosmetic preprocessing step.

This is still an exploratory result. The reported corrected AUC uses the
acquisition-day medians estimated from the full pilot. The leave-one-batch-out
control below supports batch generalization, while an independent biological
cohort is still needed before treating the result as externally validated.

Run the main check with:

```sh
cd proteotrace
python analysis_cross_batch.py \
  --matrix data/PXD001064_protein_matrix.tsv \
  --metadata data/PXD001064_sample_info.txt
```

The matrix is the pilot's target-rank-1 OpenSWATH protein-group extraction;
the metadata file is the PRIDE sample map. The script keeps one run per
declared donor and visit, selects protein groups observed in at least 99% of
runs without using donor labels, log-transforms and sample-centers them, then
subtracts the per-protein median of each acquisition day. Acquisition day and
within-day order are parsed from the run filename as a proxy because the
PXD001064 deposited files do not contain an absolute acquisition clock (see
**Timestamp check** below).

Observed output on 2026-09-19:

```text
runs_after_dedup=231 donors=116 complete_pairs=115
protein_groups_selected=273 acquisition_batches=9
raw_cross_batch_pairs: n=23483 same=93 different=23390 auc=0.7827
cross_batch_pairs: n=23483 same=93 different=23390 auc=0.9692
cross_batch_pairs: threshold_95pct_specificity=0.2115 same_pair_sensitivity=0.8602
score_vs_date_gap_same_donor: rho=-0.2611 p=0.00482 n=115
score_vs_date_gap_different_donor: rho=-0.0010 p=0.875 n=26450
date_gap_same_median=27.0 different_median=27.0 mannwhitney_p=0.403
```

The date-gap check is useful but does not replace leave-one-batch-out: same-
donor pairs are not closer in date than different-donor pairs (median 27 days
in both groups, Mann–Whitney p=0.403), while same-donor similarity still falls
with date gap (Spearman rho=-0.2611, p=0.00482).

### Leave-one-acquisition-day-out

The required batch-held-out analysis is in `analysis_leave_one_batch_out.py`:

```sh
cd proteotrace
python analysis_leave_one_batch_out.py \
  --matrix data/PXD001064_protein_matrix.tsv \
  --metadata data/PXD001064_sample_info.txt \
  --out-json data/PXD001064_leave_one_batch_out.json \
  --out-tsv data/PXD001064_leave_one_batch_out.tsv
```

Feature selection and training standardization are fitted only on the
non-held-out days. The held-out day's per-protein median is estimated from its
unlabeled runs, which matches the intended unsupervised deployment operation;
no donor label is used in preprocessing. Across the nine held-out days, the
mean fold AUC is **0.9681** and the pooled AUC after averaging the two
directions of each unique cross-day pair is **0.9705**. This supports
generalization across acquisition days, but it is not an independent biological
cohort validation and should be described as a batch-generalization control.

## MZ/DZ/unrelated control

Run:

```sh
cd proteotrace
python analysis_mz_dz_ladder.py \
  --matrix data/PXD001064_protein_matrix.tsv \
  --metadata data/PXD001064_sample_info.txt
```

The analysis applies the same preprocessing and keeps only cross-batch pairs.
The observed score ladder is:

```text
same donor median  0.3806 (n=93)
MZ twins median    0.2931 (n=80)
DZ twins median    0.1557 (n=57)
unrelated median  -0.0072 (n=23253)
```

The pairwise probabilities that the left class scores higher than the right
are 0.6273 (same donor versus MZ), 0.7224 (MZ versus DZ), 0.8276 (DZ versus
unrelated), and 0.9707 (same donor versus unrelated). The corresponding
Mann–Whitney p-values are 0.00396, 9.61e-06, 1.15e-17, and 1.66e-55.

This is a strong control-positive pattern: the score reflects relatedness in
the expected order. It also gives a required limitation for the paper. MZ
scores overlap substantially with same-donor scores, so ProteoTrace cannot be
presented as a reliable way to distinguish monozygotic twins. It is a sample
consistency check for likely swaps and duplicates, not an absolute identity
assay.

The source study is the twin plasma dataset PXD001064 ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4358658/),
[ProteomeXchange record](https://proteomecentral.proteomexchange.org/dataset/PXD001064)).

## PXD004242 external-cohort structure reconciliation

PXD004242 is a plausible independent validation cohort, but its SDRF must be
interpreted before any donor-pair metric is calculated. The reproducible check is
`analysis_pxd004242_structure.py`:

```sh
cd proteotrace
python analysis_pxd004242_structure.py \
  --sdrf data/PXD004242/PXD004242_community_annotated.sdrf.tsv \
  --manifest data/PXD004242/PXD004242_file_manifest.tsv \
  --maxquant-zip data/PXD004242/MaxQuant_Output.zip \
  --out-summary data/PXD004242/PXD004242_structure_summary.json \
  --out-samples data/PXD004242/PXD004242_biological_samples.tsv
```

The check passes with:

```text
api_files=1294 raw=1290 sdrf_rows=1290
participants=52 longitudinal=43 baseline_only=9 references=6
biological_units participant=318 longitudinal=300 baseline_only=18 references=6
participant_runs=1272 reference_runs=18
maxquant_intensity_columns=1290 label_intersection=1290
validation=PASS
```

The apparent ``58 donors`` in the SDRF is therefore not the study-donor count.
There are 52 numeric study participants and six Top14 reference-library samples
(F1--F3 and M1--M3). Of the 52 participants, 43 have longitudinal measurements
(42 with all seven timepoints and participant 44 with six), while nine have only
-8W and 0W. The 1,290 raw rows split into 1,272 participant runs and 18
reference runs. Each participant sample-timepoint has four technical runs; each
reference sample has three. This yields 318 participant biological units: 300
from the 43 longitudinal participants and 18 from the nine baseline-only
participants.

`proteinGroups.txt` contains exactly 1,290 `Intensity` columns, and every column
label maps one-to-one to an SDRF raw file. The downstream unit must be the
biological sample (`individual x sampling time`); the four technical runs are
normalized at run level for plate effects and then collapsed (for example by a
median), never counted as four independent biological observations. The six
Top14 references must be excluded from donor matching. The next analysis must
include a leave-one-plate-out check before treating an external AUC as a
final validation result. The preliminary one-run results are documented below.

The original paper reports 1,294 plasma proteomes, whereas the current PRIDE
deposit exposes 1,290 RAW files plus four non-RAW files. We keep both numbers
with their sources rather than silently treating them as the same denominator;
the SDRF/API raw-file reconciliation itself is complete.

The source is Geyer et al., *Molecular Systems Biology* 2016
([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5199119/),
[ProteomeXchange record](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD004242)).

### Layer-1 one-run validation: fixed-replicate diagnostic

The first external analysis keeps one technical run per biological sample and
uses only the 43 longitudinal participants (300 biological sample-timepoints).
It applies the same log2, per-sample centering, 99%-coverage feature selection,
and unsupervised per-batch median correction used for the PXD001064 pilot.
Because PXD001064's comparable estimate is cross-batch, the primary metric here
is cross-plate AUC; all-pair AUC is shown only as a diagnostic because many
same-donor pairs remain on the same plate.

The four fixed technical-replicate selections produced:

| selected replicate | runs | actual plates | cross-plate same-donor pairs | raw AUC | corrected AUC |
|---:|---:|---:|---:|---:|---:|
| 1 | 300 | 7 | 28 | 0.7533 | 0.7163 |
| 2 | 300 | 9 | 113 | 0.8362 | 0.9122 |
| 3 | 300 | 11 | 140 | 0.8337 | 0.9225 |
| 4 | 300 | 12 | 139 | 0.8215 | 0.8673 |

These fixed selections are a **diagnostic only**. The replicate index does not
name the same lettered plate in every plate group (for example, group 1 uses
B--E while group 2 uses A--D), and donor identity is nested in the leading-digit
plate group. Therefore the 0.7163--0.9225 range is not a statistical interval
for one estimand, and the all-pair corrected AUCs (0.9167--0.9708) must not be
compared directly with the PXD001064 cross-batch AUC 0.9692.

The extraction and analysis artifacts are:

```text
extract_pxd004242_matrix.py                 sha256 00a675c7df57723d0229558aebe562f5120235c8896966e58ce5f3e6a706bc4e
analysis_pxd004242_layer1.py                sha256 70c8153a4884b977739062196cbee49e1f29771d802b6edabe2505db3fce1340
analysis_pxd004242_leave_one_plate_out.py   sha256 eb3c8268358c74c60c2e34e33fb9bc125548cdaec6b34d6b4b8132595e4a59b2
PXD004242_protein_matrix.tsv                sha256 25a99e0f96a1841b3a405b303554412ac22158d453e0396b3e9de9680169efdf
PXD004242_layer1_fixed_summary.tsv          sha256 c9449978683a96693bb12b7ab95597894a7d6fd3a001329f45794665f4c11fd6
PXD004242_lopo_fixed_summary.tsv            sha256 52ec0f24bb352a1f45e9c0a2377b388293f49581bad50e06f253cda43abac159
analysis_pxd004242_assignment_sensitivity.py sha256 9da9eaf60b9d8345da51d9c4bb5bd45b5debe870a220c701b760c180c1b1b06
analysis_pxd004242_lopo_round.py             sha256 a6ceb801cd4e948ad1b57545d4c008229ef042c627519fbcb0e2e27cf2ef015c
PXD004242_random200_restricted.json          sha256 2f773d1e43323f4234f7b2c42a9d7e3cb4ad74dd41a70b1eea6b10a82956762f
PXD004242_random200_restricted.tsv           sha256 e91bd56bf330cad9663c2c8ea49057d7b34527c2a68c877cfe8b5dc10753f5b8
analysis_pxd004242_swap_simulation.py        sha256 7f87932c46a714d97cec21a3e2bc52c606270b9f0f074859819a796327729100
analysis_pxd004242_time_distance.py          sha256 dd60f6f40d746d70020d986df0eba1268a8984f4ddbd3ae0bc0efd5c6aa791c0
PXD004242_random200_top273.json               sha256 e9fa1dc1df4ca9659a7ca1ba0a36a38be9012d5e65b3f08f12f23be8dc943520
PXD004242_random273_panel_rotation_summary.tsv sha256 92b8498a42601cdf1184fbf928cb43e5472d4a9a0f46591ec7e92d71dc4d78a3
PXD004242_swap_within_plate.json              sha256 cca95ac6f37d647a2771c33a7576b7c2c314f927fc6b617611e0fdd56f1611f6
PXD004242_time_distance_random200.json        sha256 2fabd2dee09114622367cf9702e43e2b9ecc91aa78669ce3924a6a3cee7acf43
```

### Assignment-aware rerun

The fixed-replicate diagnostic exposed the nested design. The corrected layer-1
analysis keeps one run per sample, restricts the primary negative class to
different plates within the same leading-digit plate group, and treats the four
run choices as an assignment variable. Two checks were run:

- A timepoint-modulo-4 rotation gave corrected cross-plate AUC **0.9839**, raw
  AUC **0.8425**, sensitivity **0.9223** at 95% specificity, and **0.7760**
  at 99% specificity.
- Across **200 random assignments**, corrected AUC had median **0.9733** with
  assignment quantiles **0.9633--0.9804**; raw AUC median was **0.8463**.
  Corrected sensitivity had median **0.8768** at 95% specificity
  (quantiles 0.8408--0.9065) and **0.6967** at 99% specificity
  (quantiles 0.6483--0.7456).

These quantiles describe run-selection variability, not a statistical confidence
interval. The corrected median is close to, and slightly above, the PXD001064
0.9692 estimate; the earlier expectation of an automatic drop below 0.9692 is
not supported by this rerun. This is still an external consistency result, not
proof that the two cohorts have identical operating characteristics.

### Feature-depth check

The apparent 720-versus-273 asymmetry is a raw group count, not an effective
feature count. The adaptive 99%-coverage analysis of PXD004242 retains a median
of 147 features per assignment (range 141--153), while PXD001064 retains 272.
As an explicit cap, a fixed panel of the 273 highest-coverage PXD004242 groups
was evaluated across 200 random run assignments: corrected AUC median **0.9764**
(assignment quantiles 0.9668--0.9833), raw AUC median 0.8302, and corrected
sensitivity median 0.7154 at 99% specificity (quantiles 0.6666--0.7689).

The panel definition matters. Twenty arbitrary random 273-group panels evaluated
under the timepoint rotation had corrected AUC median **0.9527** (quantiles
0.9336--0.9679) and 99%-specificity sensitivity median 0.5264. Sparse groups
should therefore not be used to claim that feature count alone determines
transfer. The defensible result is that a high-coverage 273-group panel preserves
the transfer signal, while arbitrary low-coverage groups degrade it.

### Within-plate swap simulation

A realistic stress test exchanged the assigned donor labels of 1,000 pairs chosen
from the same actual plate, using one random run per biological sample. The
corrected adaptive panel used 148 features. Leave-two-out donor matching had
97.3% baseline top-donor accuracy and detected **99.75%** of the 2,000 swapped
queries; detection was 100% among queries whose baseline donor call was correct.
This is a simulated detection result, not a real unknown-cohort swap validation,
and it uses longitudinal reference samples that are unavailable in a truly new
cohort.

Because the 99.75% quantity is conditional on repeated donor references, the
final stress-test figure also includes a reference-count curve. The same 1,000
within-plate swaps were re-evaluated with 20 random reference subsets at each
assumed visit count. At two total visits per donor, correct recovery of the true
donor after a swap was **88.2%** (assignment quantiles 86.9--89.1%), while a
top-donor disagreement with the assigned label occurred in **99.2%** of queries.
At seven visits, correct recovery was **97.4%** and disagreement was **100%**;
the seven-visit point includes 1,958 eligible queries because one participant
has only six longitudinal timepoints. Thus disagreement is the sensitive
screening signal, whereas correct recovery is the conservative identity-aware
operating measure. Neither curve is an estimate for an unknown cohort without
longitudinal references.

The reference-count script and compact summary are
`analysis_pxd004242_swap_reference_curve.py` (sha256
`79fd875f83ea911b02e50cacd3ddd9f8d5b0b74151fae19ad5f6aad07da31494`) and
`data/PXD004242/PXD004242_swap_reference_curve_summary.json` (sha256
`e77b435b6af3179629e625b01060177a8bd59c095e1cc849693cd841f1a71916`).

### Time-distance analysis

Across the same 200 random assignments, same-donor corrected similarity was
slightly lower for pairs crossing from -8W to a post-intervention timepoint than
for pairs entirely at nonnegative timepoints: median score 0.3117 versus 0.3295,
with AUC 0.9660 versus 0.9769. The paired assignment difference in AUC had
median 0.0109 and 95% assignment quantiles -0.0031--0.0319; 177 of 200 draws
were higher for post-only pairs. At the longest 60-week gap, AUC was 0.9640
versus 0.9886 at the 4-week gap. This supports persistence of donor identity
with modest attenuation, while the per-gap estimates remain exploratory.

### Leave-one-plate-round-out transfer check

A round is the actual plate suffix (A--E), held out across all leading-digit
groups. With the timepoint rotation, pooled unique AUC was **0.9835** and the
mean fold AUC was **0.9828**. One random assignment (seed 20260919) gave pooled
AUC **0.9620** and mean fold AUC **0.9627**. This is transfer of unsupervised
plate correction; donor identity remains nested in the leading-digit groups, so
it is not donor-held-out validation.

### Leave-one-plate-out transfer check (fixed-index diagnostic)

The same four one-run selections were evaluated with feature selection and
training standardization fitted on the other plates; the held-out plate median
was estimated from its unlabeled runs. Pooled unique cross-plate AUCs were
0.7327, 0.9205, 0.9275, and 0.8790 for technical replicates 1--4, with mean
fold AUCs 0.6529, 0.8610, 0.8375, and 0.7280. Several folds in selections 1
and 2 contain only one to three held-out runs, so the pooled values are more
useful than an unweighted mean of unstable tiny folds. This is a transfer test
for unsupervised plate correction; it is not donor-held-out validation.

## Sex-check follow-up

A sex check was implemented as a separate, prespecified feasibility analysis.
It uses only the known plasma markers SHBG and PZP, with labels taken from the
external cohort's sample descriptions; the similarity score is not used to
make the sex call. Run it with:

```sh
cd proteotrace
python scripts/sex_check.py \
  --report data/PXD062541_report.tsv \
  --metadata data/PXD062541_sample_metadata.txt \
  --out-json data/PXD062541_sex_check_results.json \
  --out-tsv data/PXD062541_sex_check_values.tsv
```

The result is **not strong enough to add a production sex-QC layer**:

```text
SHBG+PZP, baseline complete: 17/28
SHBG+PZP, leave-one-participant-out: accuracy=0.7647, AUC=0.6944
SHBG+PZP, baseline model applied after intervention: 19 samples, accuracy=0.6842
PZP only, baseline complete: 26/28
PZP only, leave-one-participant-out: accuracy=0.6538, AUC=0.7143
```

The correct paper treatment is to report this as a negative feasibility result
or a future validation path, not as an independent QC feature. The PXD021677
paper reports a sex-associated plasma signal, but its deposited SDRF does not
provide usable per-sample sex labels; the PXD062541 sample-description labels
are therefore the cleaner external test used here. See the [PXD021677 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7950134/)
and [PXD062541 ProteomeXchange record](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD062541).

## Timestamp check

PXD001064 has 99 deposited raw files: 52 `.wiff`, 45 `.wiff.scan`, and only 2
`.mzXML`. The two mzXML headers were fetched by HTTP range and inspected with
`scripts/probe_timestamp_header.py`:

```sh
cd proteotrace
python scripts/probe_timestamp_header.py \
  data/PXD001064_mzxml_head_008.mzXML \
  data/PXD001064_mzxml_head_005.mzXML \
  --out data/PXD001064_timestamp_probe.tsv
```

The mzXML headers contain relative `retentionTime` values and a
`FileConverter completionTime`, but no `startTimeStamp`. The completion time is
a conversion timestamp, not the instrument acquisition time. The probe table
records the exact byte ranges, hashes, and observed fields. Consequently, the
filename date/order proxy is necessary for this validation cohort. Absolute
acquisition-time recovery remains a possible extension for repositories that
deposit mzML with an acquisition timestamp; it is not a claim of the current
PXD001064 paper.

## Inputs and hashes

```text
PXD001064_protein_matrix.tsv       44ed50af39418d3da62a96a516ba89217dbe1d1cbf2412f2db9a2f65a7878df1
PXD001064_sample_info.txt          a11c7beb614befafb5d893a9d99dc3eced29b72457ff6f16b4c72143cfd53071
PXD062541_report.tsv               88cd291cb2472635c97644191608fc73f6b1f3562fb5924d92e7e2b883832f25
PXD062541_sample_metadata.txt      5b5055660d0fdb8fd53b2523ff83b2e96dd43aa1fd9ed33627497adcc17f5b08
PXD062541.sdrf.tsv                 b9cda30fa888991dfcb9bc00a0d8130f742b686f4b7ab1019dc1319166c404a4
PXD001064_timestamp_probe.tsv      c02e26f9aa97c9bba1ade943d92ca0e7aed6d4f716159ed117ee2ed3edaea833
PXD001064_leave_one_batch_out.json 7fca178101ff6df564fdf9d575c829204b60e1ef9f2288279069be6e4bc1679b
PXD001064_leave_one_batch_out.tsv  5405a39a5add8e2717cf7ff16c512d47433d8a95355067b166310cf7305fd61d
PXD004242_community_annotated.sdrf.tsv eab1dcc720e8569ed90945b1265cf452ffe2ed5d427ffa66a26a552b23d39f69
PXD004242_file_manifest.tsv          575e0c6b8961c61009a66e195fb6864d08c7e87f6fef2d36db97685ec94034c6
PXD004242_structure_summary.json      8d0910a6810f980c878ddb2f733da77f0725e5934270d176edf64902cc5be9c9
PXD004242_biological_samples.tsv      6e83581d18e994dd1e3b1fd102b65ef51b80b4af54d4ffbb95fbbbb0b6622fde
PXD004242/MaxQuant_Output.zip         a72086382008e5059b2a465b6b1f0a651acfd09ca4cb4f58befceeff8d05e83d
```

The PXD001064 pilot's source OpenSWATH extraction was copied separately from
the deposited raw files. This workspace stores the small extracted matrix and
metadata plus the PXD062541 report used for the external feasibility check;
large vendor raw files are not duplicated here.

## Scope decision

The analysis set for the first paper is now frozen: corrected similarity,
raw versus corrected AUC, assignment-aware external validation, feature-depth
sensitivity, within-plate swap stress testing, time-distance stratification,
leave-one-plate-round-out transfer, the MZ/DZ/unrelated ladder, and transparent
limitations on sex and timestamps. Do not add a multiverse, raw-data repository
scan, or another validation cohort before writing.
