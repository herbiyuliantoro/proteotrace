# Tutorial

This tutorial uses the included PXD001064 protein-abundance matrix and sample
annotations. The example output below was obtained with the default settings.

---

## 1. Install

```sh
git clone https://github.com/herbiyuliantoro/proteotrace.git
cd proteotrace
pip install -e .
proteotrace --version
```

Python 3.10 or newer is required. pip installs numpy and pandas automatically.

---

## 2. Run the included example

The matrix and sample table for this validation cohort are included:

```sh
proteotrace check \
  --matrix  data/PXD001064_protein_matrix.tsv \
  --samples examples/PXD001064_samples.tsv \
  --out     qc/
```

The console output is:

```text
231 runs · 116 donors · 273 features
match threshold 0.3161 at 99% specificity
swap margin 0.20

  swap_suspected       7
  review              53
  unscorable           0
  no_reference        45
  ok                 126

Runs whose best match is a different donor:

                  run  claimed_donor  score_to_claimed_donor best_other_donor  score_to_best_other
yanliu_L130220_003_SW Twin pair-41-B                   0.037   Twin pair-53-B                0.473
yanliu_L130109_027_SW Twin pair-41-B                   0.037   Twin pair-41-A                0.370
yanliu_L130104_018_SW Twin pair-18-A                   0.190   Twin pair-18-B                0.452
yanliu_L130211_009_SW Twin pair-12-A                   0.167   Twin pair-44-B                0.425
yanliu_L130104_027_SW Twin pair-29-A                   0.242   Twin pair-44-B                0.480
yanliu_L130220_015_SW Twin pair-28-A                   0.200   Twin pair-28-B                0.422
yanliu_L130313_029_SW Twin pair-49-B                   0.142   Twin pair-10-B                0.343
```

---

## 3. Interpret the console output

**`273 features`** — the protein groups passing the coverage filter. Selecting
proteins measured in nearly every run provides a consistent set for comparison.
`summary.json` records that 439 proteins were excluded here.

**`threshold 0.3161`** — the 99th percentile of eligible different-donor pair
scores in this cohort, using the supplied donor annotations. The threshold is
estimated for each cohort because the score distributions can differ between
datasets. Pairwise specificity does not set the false-flag rate per run.

**`swap margin 0.20`** — the minimum difference between the best other-donor
score and the assigned-donor score for a suspected swap. The assigned-donor
score must also be below threshold and the other-donor score above threshold.
See §8 for settings.

---

## 4. Read the verdicts

| Verdict | Count | Meaning |
|:--|--:|:--|
| `ok` | 126 | Assigned-donor score reaches the match threshold |
| `review` | 53 | Assigned-donor score is below threshold, but the suspected-swap criteria are not met |
| `swap_suspected` | 7 | Assigned-donor score is below threshold; another donor reaches threshold and exceeds it by more than the margin |
| `unscorable` | 0 | Too few shared features to score the run against at least half of the other runs |
| `no_reference` | 45 | No usable repeat run from the assigned donor |

Of the 45 `no_reference` runs, 1 belongs to a donor sampled once. The other
44 belong to donors whose two visits were measured on the same acquisition day.
With multiple batches, comparisons normally use runs from different batches to
reduce the influence of shared technical effects. These donors therefore have
no eligible repeat run. When all runs are in one batch, or batch labels are
absent, the CLI uses all eligible pairs instead. If a run has no scorable
cross-batch partner at all, the report also falls back to its scorable pairs.

The unmodified validation cohorts contain 53 `review` runs out of 231 and
52 out of 300. At 99% specificity, about two in three same-donor pairs exceed
the threshold; with two visits, the assigned-donor score rests on one pair.
In PXD001064, `review` and `ok` runs have similar median missingness
(20.5% vs 18.7%) and total intensity. A `review` verdict alone does not establish
poor measurement quality or a label error. Check it alongside other QC results.

`unscorable` means that fewer than half of the possible pair scores are finite.
This can occur when too few proteins are measured in common with other runs.

---

## 5. Inspect the seven flagged runs

Check `best_other_donor` against the study metadata:

```text
yanliu_L130109_027_SW   Twin pair-41-B  →  Twin pair-41-A    ← the donor's own twin (DZ)
yanliu_L130104_018_SW   Twin pair-18-A  →  Twin pair-18-B    ← the donor's own twin (DZ)
yanliu_L130220_015_SW   Twin pair-28-A  →  Twin pair-28-B    ← the donor's own twin (MZ)

yanliu_L130220_003_SW   Twin pair-41-B  →  Twin pair-53-B    ← unrelated donor
yanliu_L130211_009_SW   Twin pair-12-A  →  Twin pair-44-B    ← unrelated donor
yanliu_L130104_027_SW   Twin pair-29-A  →  Twin pair-44-B    ← unrelated donor
yanliu_L130313_029_SW   Twin pair-49-B  →  Twin pair-10-B    ← unrelated donor
```

Three of the seven point at the donor's own twin: one identical, two
fraternal. Identical twins score almost as
high as the same person, and fraternal twins clear the threshold about one
time in five (against fewer than one in a hundred for unrelated pairs). These
three flags do not establish a sample mix-up.

The other four point at unrelated donors and need review. In the simulated-swap
evaluations, about 2–4% of runs from uninvolved donors were also flagged.
When a sample comes from someone outside the cohort, the tool still names the
best-scoring enrolled donor. That donor is a candidate to investigate; the
score does not identify the sample's source. Both runs of Twin pair-41-B are
flagged here, so neither provides a reliable reference for the other.

When reviewing a flag, check for relatives and for participants entered under
more than one study ID.

---

## 6. Output files

```text
qc/run_report.tsv    one row per run — start here
qc/pair_scores.tsv   every pair and its score — for your own plots
qc/summary.json      counts, thresholds, settings, caveats
```

`run_report.tsv` columns:

| Column | Meaning |
|:--|:--|
| `run` | run identifier, as in your matrix |
| `claimed_donor` | the donor your annotation says it is |
| `batch` | the batch you assigned |
| `n_reference_runs` | how many other runs of that donor it was compared to |
| `score_to_claimed_donor` | median score to those runs |
| `best_other_donor` | the donor that matched best instead |
| `score_to_best_other` | that donor's score |
| `n_scored_pairs` | how many other runs it could be scored against at all |
| `margin` | claimed minus best other; negative means the rival won |
| `verdict` | ok / review / swap_suspected / no_reference / unscorable |

Check `n_reference_runs`. A value of 1 means that the assigned-donor score rests
on a single pair. Additional visits provide more comparisons; detection rates
by visit count are reported in the analysis log.

---

## 7. Prepare a sample table

Three columns, tab or comma separated:

```text
run           donor      batch
Plate1_A01    subj_014   plate1
Plate1_A02    subj_014   plate1
Plate2_A01    subj_014   plate2
Plate2_A02    subj_027   plate2
```

**`run`** must match the row names of your matrix exactly.
For MaxQuant, the run names are whatever follows
`Intensity ` in `proteinGroups.txt`.

**`donor`** is the subject the sample is supposed to be from. Any string.

**`batch`** is whatever groups runs technically: plate, acquisition day, LC
column, instrument or digestion batch. Choose labels that reflect the study's
technical design. If all runs form one batch, the column may be omitted. The
CLI then uses all eligible pairs and reports the absence of batch information.
Removing batch annotation from the two development cohorts increased the
number of suspected swaps 4–10-fold; this is not an estimated error rate for
every single-batch study.

With multiple batches, the tool warns about batches containing fewer than
5 runs because their protein medians can be unstable. Batch medians are
estimated from the cohort being checked, including its evaluation runs.
Use these medians rather than estimates from another study.

For studies with multiple batches, distribute repeat samples across batches
where possible. A donor whose runs share one batch normally lacks an eligible
reference. The single-batch mode and the fallback described in §4 are exceptions.

---

## 8. Settings

```sh
proteotrace check --matrix ... --samples ... --out qc/ \
  --margin 0.15 \
  --min-coverage 0.95 \
  --specificity 0.99
```

**`--margin`** (default 0.20) controls the required advantage of the best
other-donor score. Lower margins detect more swaps but also flag more runs from
uninvolved donors. Results from `validate_tool.py`:

| margin | swaps caught | uninvolved runs flagged |
|--:|--:|--:|
| 0.05 | 74% / 99% | 14% / 11% |
| 0.20 | 70% / 98% | 3.8% / 1.8% |
| 0.25 | 67% / 96% | 2.2% / 0.8% |

(first number PXD001064, two visits per donor; second PXD004242, six to seven.)
Choose the margin according to the sensitivity and review workload needed for
the study. The default was selected on these same two development cohorts;
its performance on a new dataset can differ.

**`--rival-adjust`** (default 1, off) raises the other-donor match threshold.
In the development evaluations, it reduced sensitivity more than increasing
`--margin` for a similar reduction in flags.
Above about 10 it stops having any effect.

**`--min-coverage`** (default 0.99) decides which proteins are used: those seen
in at least this fraction of runs. A protein missing from just one run is also
kept unless the setting is 1.0. In the shallower PXD062541 cohort, reducing the
setting to 0.90 raised AUC from 0.89 to 0.94 in a sensitivity analysis.
Select the coverage rule before inspecting the verdicts to avoid choosing a
setting that favours a preferred result.

**`--specificity`** (default 0.99) sets the percentile of eligible
different-donor pair scores used as the match threshold. Increasing it raises
the threshold and usually reduces sensitivity. It does not specify the
false-flag rate per run.

---

## 9. Troubleshooting

**`only N runs appear in both files`** — your `run` column doesn't match the
matrix row names. Compare them side by side; usually it's a `.raw` extension,
a path prefix, or different capitalisation.

**`only N features are observed in 99% of runs`** — too much missingness. Try
`--min-coverage 0.95`. If that still fails, check the matrix really is runs in
rows and features in columns.

**`no positive numbers could be read`** — the file was found but none of its
values parsed as numbers. Usually the column layout or the format is wrong; try
`--format`.

**`note: values use a decimal comma`** — not an error. The export was written
with a European locale (`7,04E+12`), and the tool read the commas as decimal
points. Conversion is skipped if any value contains both `.` and `,`, or if
replacing commas does not parse at least 99% of the previously unparsed values.
Remaining unparseable values become missing values. Check the export's number
format before analysis; a mixed format does not always cause an error.

**`missing the column(s) ['donor']`** — the annotation needs `run` and `donor`
at minimum. Column names are stripped of surrounding spaces and converted to
lower case when read.

**`these runs appear more than once`** — duplicate run names. If they are
technical replicates, either give them distinct names or keep one per sample.

**`note: every run is in one batch`** — not an error, but read §7.

**`warning: N batch(es) have fewer than 5 runs`** — those runs are distorted by
their own batch median and can be flagged spuriously. Merge the batch into a
neighbouring one if that is technically defensible; otherwise discount
verdicts on those runs.

---

## 10. Limits of interpretation

- Identical twins are not reliably separable, and a twin swap can pass as
  `ok`. Fraternal twins also score high enough to appear as alternative matches.
- A single-run donor has no repeat reference. With multiple batches, donors
  whose runs share one batch normally receive `no_reference` (see §4).
- For a sample from outside the cohort, the named alternative is still an
  enrolled donor. Its identity cannot be inferred from this comparison.
- A flag does not establish the cause of an inconsistency. Possible causes
  include a swap, a mislabelled tube, or a participant entered under two IDs.
- `ok` means consistent with the assigned donor. Confirm suspected label
  errors using sample records or independent measurements.

---

Report problems through the [issue tracker](https://github.com/herbiyuliantoro/proteotrace/issues).
