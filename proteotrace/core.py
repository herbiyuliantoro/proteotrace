"""The scoring pipeline, exactly as validated on PXD001064 and PXD004242.

Nothing here uses donor labels. Labels enter only in `report.py`, to judge the
scores after the fact. Keeping that boundary in code is the reason the method
can be run on a cohort whose labels are the thing under suspicion.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# Below this many runs a batch median is mostly the runs themselves: a run alone
# in its batch is centred to zero, and two runs alone become near-perfect
# opposites (r close to -1), which can produce a spurious swap call.
MIN_BATCH_RUNS = 5


@dataclass(frozen=True)
class Normalised:
    """Feature matrix after correction, plus what was used to get there."""
    values: pd.DataFrame          # runs x features, batch-corrected and scaled
    n_features: int
    n_runs: int
    dropped_features: int
    batch_sizes: dict = field(default_factory=dict)

    @property
    def small_batches(self) -> dict:
        """Batches too small for a stable median, with their run counts."""
        return {b: n for b, n in self.batch_sizes.items() if n < MIN_BATCH_RUNS}


def normalise(matrix: pd.DataFrame, batches: pd.Series,
              min_coverage: float = 0.99) -> Normalised:
    """Select complete features, log-transform, then remove the batch offset.

    `matrix` is runs x features of positive abundances; zeros are read as
    missing. `batches` gives the batch of each run and must share its index.

    The steps are, in order: keep features observed in at least
    `min_coverage` of runs, take log2, centre each run on its own median so
    loading differences drop out, subtract the median of each feature within
    each batch, then scale each feature to unit variance so no single abundant
    protein dominates the correlation.

    A feature missing from a single run is always kept (unless `min_coverage`
    is 1). Below 100 runs a 99% rule would otherwise demand a complete
    feature, so one failed run could empty the panel. From 100 runs upward
    this changes nothing: 99% already tolerates at least one missing run.
    """
    if matrix.empty:
        raise ValueError("the feature matrix is empty")
    missing = matrix.index.difference(batches.index)
    if len(missing):
        raise ValueError(f"{len(missing)} runs have no batch assigned, "
                         f"first few: {list(missing[:5])}")
    batches = batches.loc[matrix.index]

    positive = matrix.replace(0, np.nan)
    coverage = positive.notna().mean(axis=0)
    keep = coverage >= min_coverage
    if min_coverage < 1:
        keep |= positive.isna().sum(axis=0) <= 1
    if int(keep.sum()) < 20:
        raise ValueError(
            f"only {int(keep.sum())} features are observed in {min_coverage:.0%} "
            "of runs, which is too few to score against. Lower --min-coverage, "
            "or check that the matrix is runs x features.")

    log_values = np.log2(positive.loc[:, keep])
    log_values = log_values.subtract(log_values.median(axis=1), axis=0)
    corrected = log_values - log_values.groupby(batches).transform("median")
    scaled = (corrected - corrected.mean(axis=0)) / corrected.std(axis=0).replace(0, np.nan)
    return Normalised(values=scaled, n_features=int(keep.sum()),
                      n_runs=int(matrix.shape[0]),
                      dropped_features=int((~keep).sum()),
                      batch_sizes={str(b): int(n) for b, n in batches.value_counts().items()})


def similarity(normalised: Normalised, min_overlap: int = 30) -> pd.DataFrame:
    """Pairwise correlation between runs, as a square frame."""
    return normalised.values.T.corr(min_periods=min_overlap)


def pair_table(scores: pd.DataFrame, donors: pd.Series,
               batches: pd.Series) -> pd.DataFrame:
    """Every run pair once, labelled by donor agreement and batch."""
    runs = list(scores.index)
    left_idx, right_idx = np.triu_indices(len(runs), 1)
    values = scores.to_numpy()[left_idx, right_idx]
    finite = np.isfinite(values)
    left_idx, right_idx, values = left_idx[finite], right_idx[finite], values[finite]
    left = np.asarray(runs)[left_idx]
    right = np.asarray(runs)[right_idx]
    return pd.DataFrame({
        "left": left,
        "right": right,
        "score": values,
        "same_donor": donors.loc[left].to_numpy() == donors.loc[right].to_numpy(),
        "same_batch": batches.loc[left].to_numpy() == batches.loc[right].to_numpy(),
    })


def threshold_from_pairs(pairs: pd.DataFrame, specificity: float = 0.99,
                         cross_batch_only: bool = True,
                         rivals_per_run: int = 1) -> tuple[float, float]:
    """Two thresholds from this cohort's own different-donor pairs.

    Calibrating on the cohort in hand, rather than carrying a number over from
    our validation cohorts, is deliberate: the score scale shifts with
    platform, depth and matrix, so a fixed cut-off would not transfer.

    Two thresholds are returned because they answer different questions.

    `match` is the plain quantile at `specificity`, and asks whether a run
    looks like its own donor.

    `rival` asks whether some *other* donor looks like this run, and that
    question is asked once per rival run. At 99% per pair and 300 rivals you
    would expect roughly three spurious rivals for every run, which is why a
    naive threshold buries the user in false alarms. `rival` therefore spends
    the same error budget across all the comparisons actually made. It is
    capped at the highest quantile the observed pairs can support, so it never
    extrapolates past the data; once the cap binds (about rivals_per_run > 10
    on our cohorts), raising rivals_per_run further changes nothing.

    The batch medians this score depends on must be estimated on the cohort
    being checked, including the runs under test. Medians carried over from
    other samples do not centre the new runs and inflate every score.
    """
    negatives = pairs.loc[~pairs["same_donor"]]
    if cross_batch_only and (~negatives["same_batch"]).any():
        negatives = negatives.loc[~negatives["same_batch"]]
    if negatives.empty:
        raise ValueError("no different-donor pairs, so no threshold can be set")

    values = negatives["score"].to_numpy()
    match = float(np.quantile(values, specificity))

    adjusted = 1.0 - (1.0 - specificity) / max(int(rivals_per_run), 1)
    supported = 1.0 - 20.0 / len(values)          # keep ~20 pairs in the tail
    rival = float(np.quantile(values, min(adjusted, max(supported, specificity))))
    return match, rival
