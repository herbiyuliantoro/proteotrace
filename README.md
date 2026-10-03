<div align="center">

# 🧬 ProteoTrace

**Catch swapped samples in longitudinal proteomics — from the mass spec data alone.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.12+](https://img.shields.io/badge/Python-3.12+-3776AB.svg?logo=python&logoColor=white)](requirements.txt)
[![Figures](https://github.com/herbiyuliantoro/proteotrace/actions/workflows/figures.yml/badge.svg)](https://github.com/herbiyuliantoro/proteotrace/actions/workflows/figures.yml)
[![Cohorts](https://img.shields.io/badge/validated%20on-2%20public%20cohorts-orange.svg)](#-the-data)
[![No genotype needed](https://img.shields.io/badge/genotype%20data-not%20required-brightgreen.svg)](#-why-this-exists)

</div>

---

## 🤔 Why this exists

Big clinical proteomics studies mix samples up. It happens during collection,
aliquoting, plating, injection — and most of the time nobody notices. When
someone did go looking, in the PsychENCODE BrainGVEX cohort, **19% of samples
turned out to be wrong**.

There is already a tool for this, SMAP. It works well, but it needs matching
genotype data. Most cohorts either never generated it, or can't share it.

ProteoTrace asks the same question using only what the lab already has: the
quantitative output of the run. 🎯

---

## 💡 What it does

Give it two runs and it tells you how likely they came from the same person.

The trick is not clever maths. Two unsupervised steps, neither of which looks at
donor labels, turn a weak comparison into a usable test: **scale every protein to
unit variance**, so abundant proteins stop dominating the correlation, and
**subtract each protein's median within its own batch**, so the batch offset
stops looking like shared identity.

Both are needed, and they do different jobs. Scaling carries most of the AUC;
the batch median carries the strict, high-specificity end, which is where a QC
screen operates:

```
AUC / sensitivity at 99% specificity
                 raw          scaling only   batch only     both
PXD001064   0.784 / 0.12    0.933 / 0.41   0.883 / 0.23   0.971 / 0.69
PXD004242   0.844 / 0.28    0.958 / 0.59   0.944 / 0.52   0.967 / 0.66
```

PXD004242 values are medians over 200 random one-run-per-sample assignments.

---

## 📊 Does it work?

Three independent public cohorts, three labs, three platforms. 👇 Each row uses
the same estimand: same-donor pairs against unrelated different-donor pairs, on
different batches wherever batches are known. Brackets are donor-cluster
bootstrap 95% CIs.

| Cohort | Platform | Samples | AUC | Sensitivity at 99% specificity |
|:--|:--|:--|:--:|:--:|
| 🧪 PXD001064 — twin plasma | DIA / OpenSWATH | 116 donors × 2 visits | **0.971** [0.942, 0.989] | 0.69 [0.51, 0.84] |
| 🧪 PXD004242 — weight-loss plasma | DDA / MaxQuant | 43 donors × 6–7 visits | **0.967** [0.939, 0.984] | 0.66 [0.54, 0.77] |
| 🧪 PXD062541 — epicatechin plasma ¹ | DIA / DIA-NN | 28 donors × 2 visits | **0.892** [0.792, 0.956] | 0.25 [0.07, 0.63] |

¹ The third cohort was not used to develop anything; it is the tool at its
defaults, exactly as `proteotrace check` runs it (`analysis_pxd062541.py`). No
batch information is deposited for it, so no batch correction is possible. It is
also shallower: a median of 44% of values are missing per run, and 133 proteins pass
the coverage rule. The method works there, but clearly less well. Allowing more
missing proteins raises it to 0.92–0.94 (sensitivity 0.36–0.50), a sensitivity
analysis rather than the primary result.

The operating point is the honest number to plan around, and it is uncertain.
In the two cohorts with batch information, roughly two in three genuine
same-donor pairs clear a 99%-specificity threshold, give or take 0.15. Without
batch information, and at lower depth, expect much less.

A few things worth knowing:

- 🏃 **Identity survives real biology.** The second cohort lost 12% of body
  weight over the study (Geyer et al.). Donors still matched, out to 60 weeks
  between visits: AUC 0.96 for pairs spanning the weight loss, 0.97 for pairs
  after it. This agrees with Geyer et al. (2021, *Mol Cell Proteomics*), who
  re-identified the individuals of this study from their plasma proteomes and
  pointed out the use against sample mix-ups. What ProteoTrace adds is the
  batch handling, a threshold calibrated per cohort with a measured false-alarm
  rate, the test across cohorts, and a tool.
- 🔬 **Feature quality beats feature count.** The 273 highest-coverage proteins
  give sensitivity 0.67; 273 picked at random (mostly sparse) give 0.42.
- 🔁 **Repeat visits make it much stronger.** On planted swaps, the shipped tool
  catches 60% when each donor has 2 visits and 97% with 7.
- 🧾 **Other ways to cut the data.** Including twin pairs among the negatives
  (PXD001064) gives 0.969; restricting PXD004242 negatives to one plate group
  gives 0.973. Both are in the [analysis log](docs/ANALYSIS_LOG.md).

---

## ⚠️ What it can't do

Stated up front, because you'll hit these:

- 👯 **Relatives.** MZ twin pairs score almost as high as same-donor pairs, and
  a twin swap will slip through. Dizygotic twins clear the 99% threshold 19% of
  the time, against 0.8% for unrelated pairs, so siblings and other close
  relatives in a cohort are a likely source of false rivals too.
- 📅 **Donors whose runs all share one batch.** Pairs on the same batch are
  never used, so such donors cannot be assessed. In PXD001064 that is 44 of
  the 45 unassessable runs; only one is a genuine single-run donor.
- 🎲 **Swap results are simulated.** We exchanged labels between same-batch
  runs, which is what a tube swap after plating produces. A sample from
  someone *outside* the cohort is a different error: it is usually attributed
  to an innocent enrolled donor. Nothing here is validation on a cohort where
  the truth is unknown.
- 🧪 **Batches need enough runs, and medians must come from the cohort itself.**
  Batches under 5 runs distort their runs (the tool warns). Batch medians
  carried over from other samples do not centre new runs.
- 🧫 **Plasma only, three labs.** Other matrices are untested, and the one cohort
  without batch information performs clearly worse (AUC 0.89).

Fitting the preprocessing on half the donors and scoring the other half costs
0.005–0.015 AUC, so the transductive fit is not what drives the numbers above.

---

## 🚀 Quick start

```sh
git clone https://github.com/herbiyuliantoro/proteotrace.git
cd proteotrace
pip install -e .

# try it on one of the validation cohorts
proteotrace check \
  --matrix  data/PXD001064_protein_matrix.tsv \
  --samples examples/PXD001064_samples.tsv \
  --out     qc/
```

To rebuild every figure in the paper instead:

```sh
pip install -r requirements.txt
for f in plot_fig*.py; do python "$f"; done
```

No downloads either way — both quantitative matrices are in the repo. 📈

---

## 🛠️ Using it on your own data

> 📖 **New here? Start with the [tutorial](docs/TUTORIAL.md).** It walks through
> a complete run on data already in this repository, shows the real output, and
> explains how to read a flag — including the three flags in that run that turn
> out to be identical twins rather than mix-ups.

```sh
pip install -e .

proteotrace check \
  --matrix   proteinGroups.txt \
  --samples  samples.tsv \
  --out      qc/
```

`samples.tsv` needs three columns. `batch` is whatever groups your runs
technically — plate, acquisition day, LC column:

```
run           donor      batch
Plate1_A01    subj_014   plate1
Plate1_A02    subj_014   plate1
Plate2_A01    subj_014   plate2
Plate2_A02    subj_027   plate2
```

It reads a MaxQuant `proteinGroups.txt` (or its zip), a DIA-NN report, or any
plain runs × features table, and works it out from the file. You get back a
per-run verdict:

| Verdict | What it means | What to do |
|:--|:--|:--|
| ✅ `ok` | Matches its own donor above threshold | Nothing |
| 🔍 `review` | Doesn't convincingly match its own donor, but nobody else claims it either | Expect 17–23% of runs here even in a clean cohort; look at it if other evidence points the same way |
| 🚩 `swap_suspected` | Fails its own donor, and another donor clears the threshold and beats it by the margin | Investigate before using the donor. The named rival is the best of many, not proof of the source |
| ➖ `no_reference` | Donor has no other run on a *different* batch | Out of scope; re-plate across batches if this matters |
| ⛔ `unscorable` | The run shares too few features with the others to be scored | A failed measurement; check the run |

Plus `pair_scores.tsv` and a `summary.json`.

The threshold isn't hard-coded — it's calibrated on **your** cohort's own
different-donor pairs, because the score scale shifts with platform, depth and
matrix.

### 📐 How well does it actually work?

We plant swaps we know the answer to, then check what the tool catches.
Twelve reciprocal same-batch swaps per repeat, drawn over all batches, 30
repeats, default margin 0.20. `python validate_tool.py` reproduces this:

| Cohort | Planted swaps caught | Uninvolved runs flagged | Other runs of the swapped donors flagged |
|:--|:--:|:--:|:--:|
| PXD004242 — 6–7 visits/donor, no twins | **97.6%** | 1.8% | 4.7% |
| PXD001064 — 2 visits/donor, twin cohort | **69.8%** | 3.8% | 76% |

Monte-Carlo standard errors are 2 points or less. The last column is not a
false alarm in the usual sense. When donors A and B are swapped, A's and B's
*other* runs are now compared against a wrong reference. With only two visits,
the tool can see that something is wrong among those four runs but cannot say
which one.

On the clean cohorts, with nothing planted, `swap_suspected` is 7 of 231 runs
(PXD001064) and 5 of 300 (PXD004242). `review` is 53 and 52.

The second row is the hard case on purpose: many donors have an identical twin
in the cohort, and with two visits the "does it match itself" estimate rests on
a single pair. **More visits per donor makes this much better.** Figure 6B shows
the tool catching 60% at two visits and 97% at seven.

On the third cohort (PXD062541, no batch information, two visits), 6 planted
swaps per repeat: 38% caught, 3.5% of uninvolved runs flagged; on the clean
cohort 4 of 56 runs are `swap_suspected` and 38 are `review`.
`python analysis_pxd062541.py` reproduces this.

`--margin` is the knob. Lower catches more swaps and flags more runs. 0.20 was
chosen on these same two cohorts. Tuning it on half the donors and testing on
the other half loses about 0.01–0.02 in Youden J, and a margin tuned on either
cohort lands within 0.01 of the other's optimum, so the selection is not what
drives the numbers above.

## 🗂️ What's in here

| Path | What it is |
|:--|:--|
| `proteotrace/` | 🛠️ The tool: `core`, `io`, `report`, `cli` |
| `validate_tool.py` | 📐 Plants known swaps and measures what the tool catches |
| `examples/` | 📋 Sample annotation files for the three cohorts |
| `plot_fig*.py` | 🎨 The six main figures and Figure S1, one script each |
| `proteotrace_style.py` | 🖌️ Shared style: type, colour, sizing, export |
| `figure_data.py` | 🔌 Loaders, so every figure reads the same frozen numbers |
| `analysis_*.py` | 🧮 The analyses that produced those numbers |
| `export_pxd004242_pairs.py` | ✅ Exports per-pair scores, refuses to write if they drift |
| `data/` | 💾 Matrices, sample tables, frozen result summaries — see [`data/README.md`](data/README.md) |
| `figures/` | 🖼️ PDF (vector) and PNG (600 dpi) |
| `docs/TUTORIAL.md` | 📖 Step-by-step walkthrough on real data |
| `docs/FIGURES.md` | 📝 What each figure claims, plus draft captions |
| `docs/ANALYSIS_LOG.md` | 📓 The full running record, including retired conclusions |
| `docs/AUDIT_BRIEF.md`, `docs/AUDIT_REPORT.md` | 🔎 The pre-preprint adversarial audit and what it found |
| `audit/` | 🧪 The audit's independent re-implementation (reads only raw inputs) |
| `archive/` | 🗄️ Earlier figure set, kept for provenance |

---

## 🔍 Reproducibility

We'd rather show it than claim it:

- 🤖 **CI rebuilds everything.** [`.github/workflows/figures.yml`](.github/workflows/figures.yml)
  deletes `figures/`, installs only `requirements.txt`, and rebuilds all seven on
  a clean checkout. If it can't, the build fails.
- 🔒 **The exporter checks itself.** `export_pxd004242_pairs.py` recomputes the
  raw, scaling-only and corrected AUCs. It refuses to write unless they match
  values verified by an independent implementation (`audit/audit_lib.py`) to
  within 5×10⁻⁴. A figure can't quietly drift away from the recorded analysis.
- 🔎 **Audited before submission.** Every headline number was re-derived from
  the raw inputs by code that shares nothing with the analysis scripts; see
  [`docs/AUDIT_REPORT.md`](docs/AUDIT_REPORT.md).
- 🔤 **Fonts degrade gracefully.** Arial if it's there, DejaVu Sans if not. Only
  the typeface changes.

---

## 💾 The data

Nothing new was collected. Both cohorts are public:

- 🧬 **PXD001064** — Liu et al., longitudinal twin plasma, DIA/OpenSWATH
- 🧬 **PXD004242** — Geyer et al., weight-loss plasma, DDA/MaxQuant

Two large input files are not tracked here, because they're big and you can get
them from PRIDE in one command. [`data/README.md`](data/README.md) has the URLs
and explains what's tracked and why.

---

## 📚 Citing this

A preprint is in preparation. Until it's out, cite the repository via
[`CITATION.cff`](CITATION.cff). 🙏

---

## 📄 License

MIT for the code — see [LICENSE](LICENSE). Derived tables keep the licence of
the deposition they came from.

<div align="center">

Made by **Herbi Yuliantoro** · Lehigh University 🏛️

</div>
