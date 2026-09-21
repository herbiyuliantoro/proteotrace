# ProteoTrace

**Proteomics-only sample-consistency QC for longitudinal cohorts.**

Sample swaps and donor mix-ups are common in large clinical proteomics studies
and usually go undetected: the established tool for the problem, SMAP, needs
matching genotype data, which most cohorts cannot supply or share. ProteoTrace
asks the same question from the mass-spectrometry data alone. It scores whether
two runs come from the same donor, using an unsupervised per-batch correction
and the quantitative output a laboratory already has.

Validated on two independent public cohorts on different platforms:

| Cohort | Platform | Corrected AUC | Sensitivity at 99% specificity |
|---|---|---|---|
| PXD001064 (twin plasma) | DIA / OpenSWATH | 0.9692 | 0.677 |
| PXD004242 (weight-loss plasma) | DDA / MaxQuant | 0.9733 (median of 200 assignments) | 0.697 |

Donor identity is recovered across a 12% weight-loss intervention and out to
60 weeks between visits. Batch correction is the method, not a preprocessing
step: without it the same comparison gives AUC 0.7827 and 0.8463.

**Stated limits.** Monozygotic twins are not reliably separable. All swap
results are simulations with known reference donors, not validation on a cohort
of unknown status. Donors are never held out, so these are consistency results
rather than generalisation to unseen donors.

## Quick start

```sh
git clone https://github.com/<ORG-OR-USER>/proteotrace.git
cd proteotrace
pip install -r requirements.txt
for f in plot_fig*.py; do python "$f"; done   # rebuilds every main figure
```

The two quantitative matrices are tracked, so the figures rebuild without any
download. See `data/README.md` for the two large public inputs that are not
tracked and how to fetch them.

## Reproducibility

The tracked files are enough to rebuild every main figure. This is checked two
ways, not asserted:

- `.github/workflows/figures.yml` deletes `figures/`, installs only
  `requirements.txt`, and rebuilds all six from a clean checkout on every push.
- `export_pxd004242_pairs.py` recomputes the raw and corrected AUCs of the
  assignment it exports and refuses to write unless both match the frozen
  rotation result to within 5e-4, so a figure cannot drift from the recorded
  analysis without the export failing first.

Arial is used when available and DejaVu Sans otherwise, so the figures build on
any platform; only the typeface changes.

## Paper claim currently supported

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
