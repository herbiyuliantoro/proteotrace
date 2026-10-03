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
| ⛔ `unscorable` | 0 | Shares too few features with the other runs to be scored at all |
| ➖ `no_reference` | 45 | This donor has no other run on a different batch, so nothing can be said |

`no_reference` is large here for a design reason, not because donors were
sampled once. Only 1 of these 45 runs is a genuine single-run donor. The other
44 belong to donors whose two visits were both run on the same acquisition day.
The tool never compares two runs from the same batch, because a shared batch
offset looks like shared identity, so those donors have nothing to be compared
against. If assessing every donor matters, spread each donor's runs over
different batches.

`review` is not a soft accusation, and it is common: 53 of 231 here, 52 of 300
in the other validation cohort, with nothing planted. It mostly reflects the
ceiling of the score itself. At 99% specificity, only about two in three genuine
same-donor pairs clear the threshold, and with two visits the judgement rests on
one pair. Review runs here are barely different in quality from `ok` runs
(median missingness 20.5% vs 18.7%, same total intensity). So treat `review`
as "not confirmed", not as "bad run", and look harder only if something else
points the same way.

`unscorable` is the verdict for a run that failed as a measurement. If most of
its values are missing, it cannot be scored against at least half of the other
runs.

---

## 5. Now look harder at those 7 flags 🧐

This is the part that matters, and it is why the tool prints the rival rather
than just a verdict. Look at who the rival is:

```text
yanliu_L130109_027_SW   Twin pair-41-B  →  Twin pair-41-A    ← the donor's own twin (DZ)
yanliu_L130104_018_SW   Twin pair-18-A  →  Twin pair-18-B    ← the donor's own twin (DZ)
yanliu_L130220_015_SW   Twin pair-28-A  →  Twin pair-28-B    ← the donor's own twin (MZ)

yanliu_L130220_003_SW   Twin pair-41-B  →  Twin pair-53-B    ← unrelated donor
yanliu_L130211_009_SW   Twin pair-12-A  →  Twin pair-44-B    ← unrelated donor
yanliu_L130104_027_SW   Twin pair-29-A  →  Twin pair-44-B    ← unrelated donor
yanliu_L130313_029_SW   Twin pair-49-B  →  Twin pair-10-B    ← unrelated donor
```

Three of the seven point at the donor's **own twin**: one identical, two
fraternal. That is expected in a twin cohort. Identical twins score almost as
high as the same person, and fraternal twins clear the threshold about one
time in five (against fewer than one in a hundred for unrelated pairs). These
three are not evidence of a mix-up.

The other four point at **unrelated** donors. Those are worth a look, but keep
the base rate in mind. On a clean cohort, the tool flags about 2–4% of runs
this way by chance. And when a sample comes from someone outside the cohort,
the tool still names *some* enrolled donor as the best rival. A flag says "this
run does not look like its claimed donor"; the named rival is a lead, not an
identification. Two flags on the same donor (both runs of Twin pair-41-B here)
are more telling than one.

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
| `n_scored_pairs` | how many other runs it could be scored against at all |
| `margin` | claimed minus best other; negative means the rival won |
| `verdict` | ok / review / swap_suspected / no_reference / unscorable |

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
say so. Expect many more false swap flags, though: without the batch
correction our two validation cohorts gave 4–10 times as many.

Two more things about batches. Each needs enough runs: below 5, the batch
median is mostly the runs themselves, and the tool warns. And the tool
deliberately computes each batch's median from *your* runs, including the
ones it is checking. Don't try to save time by reusing medians from another
study; they do not centre new runs and every score inflates.

Spread each donor's runs over different batches if you can. Pairs on the same
batch are never used, so a donor whose runs share one batch gets
`no_reference`.

---

## 8. Turning the knobs 🎛️

```sh
proteotrace check --matrix ... --samples ... --out qc/ \
  --margin 0.15 \
  --min-coverage 0.95 \
  --specificity 0.99
```

**`--margin`** (default 0.20) is the main one. It trades catching swaps against
disturbing runs that have nothing to do with the swap. Measured on our two
cohorts with `validate_tool.py`:

| margin | swaps caught | uninvolved runs flagged |
|--:|--:|--:|
| 0.05 | 74% / 99% | 14% / 11% |
| 0.20 | 70% / 98% | 3.8% / 1.8% |
| 0.25 | 67% / 96% | 2.2% / 0.8% |

(first number PXD001064, two visits per donor; second PXD004242, six to seven.)
Lower it if a missed swap is worse than a wasted hour; raise it if you are
drowning in flags. The default was picked on these same two cohorts, so treat
it as a reasonable start rather than as tuned for your data.

**`--rival-adjust`** (default 1, off) raises the bar for the rival instead. It
costs far more sensitivity than raising `--margin` for the same drop in flags.
Above about 10 it stops having any effect.

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

**`warning: N batch(es) have fewer than 5 runs`** — those runs are distorted by
their own batch median and can be flagged spuriously. Merge the batch into a
neighbouring one if that is technically defensible; otherwise discount
verdicts on those runs.

---

## 10. What it cannot tell you 🚧

- 👯 **Relatives.** Identical twins are not separable, and a twin swap can pass
  as `ok`. Fraternal twins (and so, most likely, siblings) also score high
  enough to appear as rivals.
- 🧍 **Donors with no cross-batch pair.** Single-run donors, and donors whose
  runs all share one batch, get no verdict.
- 🕵️ **Who a stray sample belongs to.** If a sample comes from someone outside
  the cohort, the tool still names the closest enrolled donor.
- ❓ **Why** a sample is wrong. It flags an inconsistency; the cause could be a
  swap, a mislabelled tube, a duplicated subject ID, or a sample from the same
  person entered twice under different IDs.
- ✋ **Proof.** `ok` means consistent with the claimed donor, not identity
  confirmed. Treat flags as a worklist, not a verdict.

---

Something unclear or broken? Open an issue:
<https://github.com/herbiyuliantoro/proteotrace/issues> 🙏
