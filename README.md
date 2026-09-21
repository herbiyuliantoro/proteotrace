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

The trick is not clever maths. It's that **batch effects drown the signal, and
removing them is the whole method**. Subtract each protein's median within its
own batch — no donor labels, nothing supervised — and a comparison that was
barely better than a coin flip becomes a usable test.

```
                raw   →   batch-corrected
PXD001064     0.7827  →      0.9692        (AUC)
PXD004242     0.8463  →      0.9733
```

---

## 📊 Does it work?

Two independent public cohorts, two different platforms, same answer. 👇

| Cohort | Platform | Samples | AUC | Sensitivity at 99% specificity |
|:--|:--|:--|:--:|:--:|
| 🧪 PXD001064 — twin plasma | DIA / OpenSWATH | 116 donors | **0.9692** | 0.677 |
| 🧪 PXD004242 — weight-loss plasma | DDA / MaxQuant | 43 donors × 7 visits | **0.9733** | 0.697 |

A few things worth knowing:

- 🏃 **Identity survives real biology.** The second cohort lost 12% of body
  weight over the study. Donors still matched, out to 60 weeks between visits.
- 🔬 **Feature quality beats feature count.** 273 well-measured proteins give
  sensitivity 0.715. The same number picked at random gives 0.526.
- 🔁 **Repeat visits make it much stronger.** With 2 visits per donor it names
  the right donor 88% of the time. With 7, it's 97%.

---

## ⚠️ What it can't do

Stated up front, because you'll hit these:

- 👯 **Identical twins.** MZ twin pairs score almost as high as same-donor
  pairs. A twin swap will slip through.
- 🎲 **Swap results are simulated.** We permuted labels on known samples. This
  is not validation on a cohort where the truth is unknown.
- 🚪 **Donors are never held out.** These are consistency results, not proof
  that the score generalises to donors the method has never seen.
- 🧫 **Plasma only, two labs.** Other matrices and other labs are untested.

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
| 🔍 `review` | Doesn't convincingly match its own donor, but nobody else claims it either | Check the run quality first |
| 🚩 `swap_suspected` | Fails its own donor, and another donor matches it better | Investigate before using the donor |
| ➖ `no_reference` | Donor has only this one run, so nothing can be said | Out of scope |

Plus `pair_scores.tsv` and a `summary.json`.

The threshold isn't hard-coded — it's calibrated on **your** cohort's own
different-donor pairs, because the score scale shifts with platform, depth and
matrix.

### 📐 How well does it actually work?

We plant swaps we know the answer to, then check what the tool catches.
`python validate_tool.py` reproduces this:

| Cohort | Catches planted swaps | Flags clean runs |
|:--|:--:|:--:|
| PXD004242 — 7 visits/donor, no twins | **93.3%** | 4.9% |
| PXD001064 — 2 visits/donor, twin cohort | **80.6%** | 12.2% |

The second row is the hard case on purpose: half the cohort has an identical
twin in it, and each donor has only two visits, so the "does it match itself"
estimate rests on a single pair. **More visits per donor makes this much
better** — that's the same dependence Figure 6B shows.

`--margin` is the knob. Lower catches more swaps and flags more clean runs;
0.20 is the default because it sat at the best point on both cohorts. The
numbers above are for the default.

## 🗂️ What's in here

| Path | What it is |
|:--|:--|
| `proteotrace/` | 🛠️ The tool: `core`, `io`, `report`, `cli` |
| `validate_tool.py` | 📐 Plants known swaps and measures what the tool catches |
| `examples/` | 📋 Sample annotation files for both cohorts |
| `plot_fig*.py` | 🎨 The six main figures, one script each |
| `proteotrace_style.py` | 🖌️ Shared style: type, colour, sizing, export |
| `figure_data.py` | 🔌 Loaders, so every figure reads the same frozen numbers |
| `analysis_*.py` | 🧮 The analyses that produced those numbers |
| `export_pxd004242_pairs.py` | ✅ Exports per-pair scores, refuses to write if they drift |
| `data/` | 💾 Matrices, sample tables, frozen result summaries — see [`data/README.md`](data/README.md) |
| `figures/` | 🖼️ PDF (vector) and PNG (600 dpi) |
| `docs/FIGURES.md` | 📝 What each figure claims, plus draft captions |
| `docs/ANALYSIS_LOG.md` | 📓 The full running record, including retired conclusions |
| `archive/` | 🗄️ Earlier figure set, kept for provenance |

---

## 🔍 Reproducibility

We'd rather show it than claim it:

- 🤖 **CI rebuilds everything.** [`.github/workflows/figures.yml`](.github/workflows/figures.yml)
  deletes `figures/`, installs only `requirements.txt`, and rebuilds all six on
  a clean checkout. If it can't, the build fails.
- 🔒 **The exporter checks itself.** `export_pxd004242_pairs.py` recomputes both
  AUCs and refuses to write unless they match the frozen result to within
  5×10⁻⁴. A figure can't quietly drift away from the recorded analysis.
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
