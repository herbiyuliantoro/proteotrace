# ProteoTrace

ProteoTrace compares protein-abundance profiles across longitudinal plasma
proteomics runs. It flags runs whose profiles differ from those of their
assigned donors. A flag needs review and does not identify a sample's source.

[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](pyproject.toml)
[![MIT license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Figure build](https://github.com/herbiyuliantoro/proteotrace/actions/workflows/figures.yml/badge.svg)](https://github.com/herbiyuliantoro/proteotrace/actions/workflows/figures.yml)

## What it does

The software takes an abundance matrix with runs in rows and proteins in
columns, plus a sample table with run, donor, and (when available) batch labels.
It log2-transforms positive abundances, centres each run on its median,
subtracts each protein's median within each batch, scales each protein to unit
variance, and calculates pairwise Pearson correlations. Donor labels do not
enter the score. They are used to set a cohort-specific threshold and compare
runs with their assigned donors.

Protein scaling accounts for most of the AUC gain. Batch centring improves
sensitivity at high specificity. The table shows AUC and sensitivity at 99%
pairwise specificity for the two cohorts with batch information. PXD004242
values are medians over 200 random run assignments.

| Cohort | Raw | Scaling only | Batch centring only | Both steps |
|:--|--:|--:|--:|--:|
| PXD001064 | 0.784 / 0.12 | 0.933 / 0.41 | 0.883 / 0.23 | **0.971 / 0.69** |
| PXD004242 | 0.844 / 0.28 | 0.958 / 0.59 | 0.944 / 0.52 | **0.967 / 0.66** |

Each cell gives AUC / sensitivity at 99% specificity. Full results and
confidence intervals are in the [analysis log](docs/ANALYSIS_LOG.md).

## Results across cohorts

The primary comparison is between same-donor pairs and pairs from unrelated
donors. For cohorts with batch labels, both pair types span different batches.
PXD062541 has no deposited batch labels, so all eligible pairs were used and no
batch centring was possible. Intervals are 95% donor-cluster bootstrap
confidence intervals.

| Cohort | Data | Runs | AUC | Sensitivity at 99% specificity |
|:--|:--|--:|--:|--:|
| PXD001064 | Twin plasma, DIA/OpenSWATH | 231 | 0.971 [0.942, 0.989] | 0.69 [0.51, 0.84] |
| PXD004242 | Weight-loss plasma, DDA/MaxQuant | 300 | 0.967 [0.939, 0.984] | 0.66 [0.54, 0.77] |
| PXD062541 | Epicatechin plasma, DIA-NN | 56 | 0.892 [0.792, 0.956] | 0.25 [0.07, 0.63] |

The third cohort was evaluated with the default settings and was not used to
develop the method. It had no batch information, a median of 44% missing values
per run, and 133 proteins passing the default coverage filter. Its results were
lower than those for the first two cohorts. The [analysis log](docs/ANALYSIS_LOG.md)
reports the coverage sensitivity analysis.

Geyer et al. (2021, *Molecular & Cellular Proteomics*) report an average 12%
weight loss in PXD004242. Same-donor AUC was 0.964 for the longest, 60-week gap
and 0.989 for the 4-week gap (analysis in this repository).

We simulated reciprocal swaps by exchanging the labels of two same-batch runs.
The tool flagged 69.8% of swapped runs in
PXD001064 and 97.6% in PXD004242. It also flagged 3.8% and 1.8% of runs from
unaffected donors, respectively. In the two-visit twin cohort, it also flagged
76% of the affected donors' other runs. After a swap, those runs no longer have
a reliable same-donor reference. The [analysis log](docs/ANALYSIS_LOG.md) gives
the simulation design and counts.

## Limits

- The method cannot distinguish monozygotic twins. Dizygotic twins also match
  above the 99% threshold more often than unrelated donors (19% versus 0.8% in
  PXD001064).
- With multiple batches, comparisons normally use runs from different batches.
  A donor whose runs all fall in one batch usually lacks an eligible reference.
  With one batch or no batch labels, the CLI uses all eligible pairs. The report
  also falls back to all scorable pairs if a run has no scorable cross-batch
  partner. These exceptions are described in the [tutorial](docs/TUTORIAL.md).
  Batch medians come from the cohort being analysed; batches with fewer than
  five runs may have unstable medians.
- Swap detection was evaluated by exchanging donor labels between runs from
  the same batch. The study did not include confirmed real sample swaps or
  samples from people outside the enrolled cohorts.
- The three evaluated cohorts are plasma studies. Performance on other sample
  types has not been established. Performance was lower in the cohort without
  batch labels and with shallower measurements.
- Thresholds depend on the cohort and its donor annotations. A suspected swap
  is a reason to investigate alongside other evidence, not proof of identity
  or sample origin.

## Install and run

ProteoTrace requires Python 3.10 or newer. Install it from a local checkout:

```sh
git clone https://github.com/herbiyuliantoro/proteotrace.git
cd proteotrace
python -m pip install -e .
```

Run the included PXD001064 example:

```sh
proteotrace check \
  --matrix data/PXD001064_protein_matrix.tsv \
  --samples examples/PXD001064_samples.tsv \
  --out qc/
```

The command writes `run_report.tsv`, `pair_scores.tsv`, and `summary.json` to
`qc/`. Read the [tutorial](docs/TUTORIAL.md) for the example output, input
formats, sample-table requirements, verdicts, and threshold settings.

To regenerate the figures, install the analysis dependencies and run the plot
scripts:

```sh
python -m pip install -r requirements.txt
for file in plot_fig*.py; do python "$file"; done
```

The included example and figure scripts run from the files tracked in this
repository. Reproducing the underlying analyses may require larger source
exports from PRIDE; download instructions and file provenance are in
[`data/README.md`](data/README.md).

## Repository contents

| Path | Contents |
|:--|:--|
| `proteotrace/` | Python package and command-line interface |
| `examples/` | Sample annotations for the three cohorts |
| `data/` | Quantitative matrices, sample metadata, and analysis summaries |
| `plot_fig*.py`, `figure_data.py`, `proteotrace_style.py` | Figure-generation code |
| `analysis_*.py`, `validate_tool.py` | Cohort analyses and simulated-swap evaluation |
| `audit/` | Independent reimplementation used for the preprint audit |
| `figures/` | Main and supplementary figures in PDF and PNG formats |
| `docs/TUTORIAL.md` | User guide |
| `docs/FIGURES.md` | Figure descriptions and captions |
| `docs/ANALYSIS_LOG.md` | Analysis history and corrections |
| `docs/AUDIT_SCOPE.md` | Scope of the independent audit |
| `docs/AUDIT_REPORT.md` | Independent audit findings |

## Data and citation

The cohort data are public under PRIDE accessions PXD001064, PXD004242, and
PXD062541. No new data were collected for this project. Derived tables retain
the terms of their source depositions; the software is licensed under MIT. See
[`data/README.md`](data/README.md) for the sources and tracked files.

To cite the software, use the repository metadata in [`CITATION.cff`](CITATION.cff).
