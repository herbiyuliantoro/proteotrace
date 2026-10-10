"""Summarise pairwise similarities and assign per-run QC verdicts."""
from __future__ import annotations

import numpy as np
import pandas as pd

VERDICTS = ("ok", "review", "swap_suspected", "no_reference", "unscorable")

# A run that cannot be scored against at least this share of the other runs
# (too few shared features, see --min-overlap) is marked unscorable.
MIN_SCORED_FRACTION = 0.5


def run_report(scores: pd.DataFrame, donors: pd.Series, batches: pd.Series,
               threshold: float, margin: float = 0.05,
               cross_batch_only: bool = True,
               rival_threshold: float | None = None) -> pd.DataFrame:
    """Compare each run with its assigned donor and best other-donor match.

    The assigned-donor score is the median over usable repeat runs. The
    other-donor score is the maximum over usable runs from other donors.
    With cross_batch_only=True, use cross-batch comparisons unless none are
    scorable for this run; then fall back to all scorable pairs.

    Verdicts:
      ok              assigned-donor score is at least threshold
      review          assigned-donor score is below threshold, but the
                      suspected-swap criteria are not met
      swap_suspected  assigned-donor score is below threshold; the best
                      other-donor score is at least rival_threshold and
                      exceeds the assigned-donor score by more than margin
      no_reference    no usable repeat run from the assigned donor
      unscorable      fewer than MIN_SCORED_FRACTION of possible pair scores
                      are finite, before applying the batch restriction

    A suspected swap does not identify the sample's source. Related donors
    can produce high scores, and a sample from outside the cohort still has
    a best-scoring enrolled donor. A review verdict alone does not establish
    poor measurement quality.
    """
    if rival_threshold is None:
        rival_threshold = threshold
    runs = list(scores.index)
    values = scores.to_numpy(dtype=float)
    donor_of = donors.loc[runs].to_numpy()
    batch_of = batches.loc[runs].to_numpy()
    position = {run: i for i, run in enumerate(runs)}

    rows = []
    for run in runs:
        i = position[run]
        row = values[i].copy()
        row[i] = np.nan
        n_scored = int(np.isfinite(row).sum())
        scorable = n_scored >= MIN_SCORED_FRACTION * (len(runs) - 1)
        usable = np.isfinite(row)
        if cross_batch_only:
            usable &= batch_of != batch_of[i]
            if not usable.any():           # single-batch donor: fall back
                usable = np.isfinite(row)
                usable[i] = False

        own = usable & (donor_of == donor_of[i])
        other = usable & (donor_of != donor_of[i])

        own_score = float(np.median(row[own])) if own.any() else np.nan
        if other.any():
            best = int(np.nanargmax(np.where(other, row, -np.inf)))
            best_score, best_donor = float(row[best]), donor_of[best]
        else:
            best_score, best_donor = np.nan, None

        claimed_ok = np.isfinite(own_score) and own_score >= threshold
        rival_convincing = (np.isfinite(best_score) and best_score >= rival_threshold
                            and best_score > own_score + margin)
        if not scorable:
            verdict = "unscorable"
        elif not own.any():
            verdict = "no_reference"
        elif claimed_ok:
            verdict = "ok"
        elif rival_convincing:
            verdict = "swap_suspected"
        else:
            verdict = "review"

        rows.append({
            "run": run,
            "claimed_donor": donor_of[i],
            "batch": batch_of[i],
            "n_reference_runs": int(own.sum()),
            "n_scored_pairs": n_scored,
            "score_to_claimed_donor": own_score,
            "best_other_donor": best_donor,
            "score_to_best_other": best_score,
            "margin": (own_score - best_score) if np.isfinite(best_score) else np.nan,
            "verdict": verdict,
        })

    report = pd.DataFrame(rows)
    order = {v: i for i, v in enumerate(("swap_suspected", "review", "unscorable",
                                         "no_reference", "ok"))}
    return report.sort_values(["verdict", "margin"],
                              key=lambda c: c.map(order) if c.name == "verdict" else c
                              ).reset_index(drop=True)


def summarise(report: pd.DataFrame, threshold: float, n_features: int,
              specificity: float, rival_threshold: float | None = None) -> dict:
    counts = report["verdict"].value_counts().to_dict()
    return {
        "runs": int(len(report)),
        "donors": int(report["claimed_donor"].nunique()),
        "features_used": int(n_features),
        "specificity_target": float(specificity),
        "score_threshold": float(threshold),
        "rival_threshold": float(rival_threshold if rival_threshold is not None else threshold),
        "verdicts": {v: int(counts.get(v, 0)) for v in VERDICTS},
        "caveats": [
            "Identical twins are not reliably separable; a twin swap can pass as ok.",
            "Other relatives also score high: in PXD001064, 19% of dizygotic-twin pairs "
            "clear the 99% threshold against 0.8% of unrelated pairs.",
            "The donor named in swap_suspected is the best-scoring rival, not an "
            "identification of the sample's source.",
            "The threshold is calibrated on this cohort's own different-donor pairs.",
            "A verdict of ok means consistent with the claimed donor, not proof of identity.",
        ],
    }
