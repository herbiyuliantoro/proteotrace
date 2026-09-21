"""Turn pair scores into a per-run verdict a person can act on."""
from __future__ import annotations

import numpy as np
import pandas as pd

VERDICTS = ("ok", "review", "swap_suspected", "no_reference")


def run_report(scores: pd.DataFrame, donors: pd.Series, batches: pd.Series,
               threshold: float, margin: float = 0.05,
               cross_batch_only: bool = True,
               rival_threshold: float | None = None) -> pd.DataFrame:
    """One row per run: how well it matches its own donor, and its best rival.

    For each run we take the median score to the other runs of its *claimed*
    donor, and the best score to any run of any other donor. Comparing those
    two is what separates "this sample is weak" from "this sample belongs to
    somebody else".

    Calling a swap takes three things at once, not one. A rival scoring
    slightly higher is normal: with only two visits per donor the own-donor
    estimate rests on a single pair, and identical twins sit almost on top of
    their sibling. So a swap is only claimed when the run fails to match its
    own donor, *and* some other donor matches it convincingly, *and* the gap
    between them is bigger than `margin`. Anything weaker is reported as
    something to look at, not as a swap.

    Verdicts:
      ok              matches its claimed donor above threshold
      review          does not match its claimed donor above threshold, but no
                      other donor matches it convincingly either
      swap_suspected  fails its claimed donor, and another donor clears the
                      threshold and beats it by more than `margin`
      no_reference    the claimed donor has no other run to compare against,
                      so nothing can be said
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
        if not own.any():
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
            "score_to_claimed_donor": own_score,
            "best_other_donor": best_donor,
            "score_to_best_other": best_score,
            "margin": (own_score - best_score) if np.isfinite(best_score) else np.nan,
            "verdict": verdict,
        })

    report = pd.DataFrame(rows)
    order = {v: i for i, v in enumerate(("swap_suspected", "review", "no_reference", "ok"))}
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
            "The threshold is calibrated on this cohort's own different-donor pairs.",
            "A verdict of ok means consistent with the claimed donor, not proof of identity.",
        ],
    }
