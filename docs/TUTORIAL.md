# 📖 Tutorial

Everything below is a real run on data that ships with this repository, so you
can follow along and get the same numbers.

---

## 1. Install 📦

```sh
git clone https://github.com/herbiyuliantoro/proteotrace.git
cd proteotrace
pip install -e .
proteotrace --version
```

You need Python 3.10 or newer, and numpy and pandas, which pip pulls in.

---

## 2. Run it on data that is already here ▶️

No download needed — one of the validation cohorts is in the repository:

```sh
proteotrace check \
  --matrix  data/PXD001064_protein_matrix.tsv \
  --samples examples/PXD001064_samples.tsv \
  --out     qc/
```

It takes about a minute. Here is what it prints:

```text
231 runs · 116 donors · 273 features
match threshold 0.3161 at 99% specificity
swap margin 0.20

  swap_suspected       7
  review              53
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

## 3. What the header is telling you 🔎

**`273 features`** — the proteins measured in nearly every run. Everything else
is dropped, because a protein missing from half the runs cannot compare them.
`summary.json` records that 439 proteins were dropped here.

**`threshold 0.3161`** — nothing is hard-coded. The tool takes the pairs in
*your* cohort that are known to be different donors, and puts the line where
only 1% of them fall above it. On a different platform the number will come
out different, and that is the point: score scales do not transfer.

**`swap margin 0.20`** — how much better another donor has to look before the
tool will say "this belongs to somebody else". See §6.

---

## 4. Reading the verdicts 🚦

| Verdict | Count | Meaning |
|:--|--:|:--|
| ✅ `ok` | 126 | Matches its own donor above threshold |
| 🔍 `review` | 53 | Doesn't convincingly match its own donor, but nobody else claims it |
| 🚩 `swap_suspected` | 7 | Fails its own donor, and another donor matches it better |
| ➖ `no_reference` | 45 | This donor has no other run, so nothing can be said |

`no_reference` is large here because many donors were sampled only once. The
tool needs at least two runs from a donor to say anything, so single-run donors
are simply out of scope. That is a property of your study design, not a fault.

`review` is not a soft accusation. It usually means the run itself is weak —
low coverage, a poor injection, degraded sample. Check the run quality before
suspecting the label.

---

## 5. Now look harder at those 7 flags 🧐

This is the part that matters, and it is why the tool prints the rival rather
than just a verdict. Look at who the rival is:

```text
yanliu_L130109_027_SW   Twin pair-41-B  →  Twin pair-41-A    ← the donor's own twin
yanliu_L130104_018_SW   Twin pair-18-A  →  Twin pair-18-B    ← the donor's own twin
yanliu_L130220_015_SW   Twin pair-28-A  →  Twin pair-28-B    ← the donor's own twin

yanliu_L130220_003_SW   Twin pair-41-B  →  Twin pair-53-B    ← unrelated donor
yanliu_L130211_009_SW   Twin pair-12-A  →  Twin pair-44-B    ← unrelated donor
yanliu_L130104_027_SW   Twin pair-29-A  →  Twin pair-44-B    ← unrelated donor
yanliu_L130313_029_SW   Twin pair-49-B  →  Twin pair-10-B    ← unrelated donor
```

Three of the seven point at the donor's **own twin**. In a twin cohort that is
expected: identical twins score almost as high as the same person, so the tool
cannot tell them apart, and it says so. Those three are a coin flip, not
evidence of a mix-up.

The other four point at **unrelated** donors. Those are the ones worth
investigating: a plasma sample matching a stranger better than itself is not
something noise usually produces.

So the practical rule: 🔑 **read the rival, not just the flag.** A rival who is
a relative, a twin, or the same person under a second study ID is a different
situation from a rival who is a stranger.

---

## 6. The three output files 📂

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
| `margin` | claimed minus best other; negative means the rival won |
| `verdict` | ok / review / swap_suspected / no_reference |

⚠️ Watch `n_reference_runs`. If it is 1, the whole judgement rests on a single
pair and is noisy. Two visits per donor is the hardest case; more visits makes
the tool much more reliable.

---

## 7. Preparing your own `samples.tsv` ✍️

Three columns, tab or comma separated:

```text
run           donor      batch
Plate1_A01    subj_014   plate1
Plate1_A02    subj_014   plate1
Plate2_A01    subj_014   plate2
Plate2_A02    subj_027   plate2
```

**`run`** must match the row names of your matrix *exactly*. This is the most
common thing to get wrong. For MaxQuant, the run names are whatever follows
`Intensity ` in `proteinGroups.txt`.

**`donor`** is the subject the sample is supposed to be from. Any string.

**`batch`** is whatever groups runs technically: plate, acquisition day, LC
column, instrument, digestion batch. Pick the one that changes most between
runs. If you genuinely have one batch, leave the column out and the tool will
say so — but scores will be optimistic, because part of what it usually
removes is not there to remove.

---

## 8. Turning the knobs 🎛️

```sh
proteotrace check --matrix ... --samples ... --out qc/ \
  --margin 0.15 \
  --min-coverage 0.95 \
  --specificity 0.99
```

**`--margin`** (default 0.20) is the main one. It trades catching swaps against
disturbing clean runs. Measured on our two cohorts:

| margin | swaps caught | clean runs flagged |
|--:|--:|--:|
| 0.05 | 86% / 97% | 22% / 16% |
| 0.20 | 81% / 93% | 12% / 5% |
| 0.25 | 81% / 92% | 11% / 3% |

(first number PXD001064, second PXD004242; the 0.20 row is what
`validate_tool.py` reports for the shipped default). Lower it if a missed swap is worse
than a wasted hour; raise it if you are drowning in flags.

**`--min-coverage`** (default 0.99) is worth lowering to 0.95 if you have deep
data with a lot of missingness and end up with too few features.

**`--specificity`** (default 0.99) moves the match threshold. Leave it alone
unless you know why you are changing it.

---

## 9. When it complains 🧯

**`only N runs appear in both files`** — your `run` column doesn't match the
matrix row names. Compare them side by side; usually it's a `.raw` extension,
a path prefix, or different capitalisation.

**`only N features are observed in 99% of runs`** — too much missingness. Try
`--min-coverage 0.95`. If that still fails, check the matrix really is runs in
rows and features in columns.

**`missing the column(s) ['donor']`** — the annotation needs `run` and `donor`
at minimum, lower-case.

**`these runs appear more than once`** — duplicate run names. If they are
technical replicates, either give them distinct names or keep one per sample.

**`note: every run is in one batch`** — not an error, but read §7.

---

## 10. What it cannot tell you 🚧

- 👯 **Identical twins.** They are not separable. A twin swap can pass as `ok`.
- 🧍 **Single-run donors.** No reference, no verdict.
- ❓ **Why** a sample is wrong. It flags an inconsistency; the cause could be a
  swap, a mislabelled tube, a duplicated subject ID, or a sample from the same
  person entered twice under different IDs.
- ✋ **Proof.** `ok` means consistent with the claimed donor, not identity
  confirmed. Treat flags as a worklist, not a verdict.

---

Something unclear or broken? Open an issue:
<https://github.com/herbiyuliantoro/proteotrace/issues> 🙏
