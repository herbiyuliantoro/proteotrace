# TASK: adversarial audit of ProteoTrace

You are auditing a proteomics QC method and tool before its author submits a
preprint. Your job is to find what is wrong with it, not to confirm it works.
The author wants to be told now, by you, rather than later by a reviewer.

Repository: /mnt/c/ProteoTrace  (private; branch main)
Run Python with: sh /mnt/c/Metabolomics_project/scripts/py <script.py>
Read docs/ANALYSIS_LOG.md for the full history, docs/FIGURES.md for the
figures, docs/TUTORIAL.md for the tool, README.md for the headline claims.

The repository is private, so the audit must run locally with filesystem
access on the author's machine.

## THE CLAIMS YOU ARE AUDITING

| # | Claim | Number |
|---|---|---|
| 1 | Batch-corrected similarity separates same-donor from different-donor pairs | AUC 0.9692 (PXD001064), 0.9733 (PXD004242) |
| 2 | Batch correction is the method, not preprocessing | 0.7827 → 0.9692 and 0.8463 → 0.9733 |
| 3 | Score tracks relatedness | same donor 0.381 > MZ 0.293 > DZ 0.156 > unrelated −0.007 |
| 4 | Operating point | sensitivity 0.677 / 0.697 at 99% specificity |
| 5 | Identity survives a 12% weight-loss intervention over 60 weeks | gap-wise AUC 0.96–0.99 |
| 6 | Coverage quality beats panel size | top-273 0.976/0.715 vs random-273 0.953/0.526 |
| 7 | The shipped tool catches planted swaps | 93.3% at 4.9% false alarms; 80.6% at 12.2% |

## ATTACK THESE FIRST — ranked by how much I expect them to hurt

1. **The tool's default margin (0.20) was tuned on the same two cohorts whose
   performance is then reported.** Claim 7 is therefore selected-on, and the
   true out-of-sample numbers are worse by an unknown amount. Quantify it.
   Nested resampling, a held-out cohort split, anything defensible. If you can
   show how much optimism this buys, that number must go in the paper.

2. **Nothing is ever held out.** Feature selection (99% coverage), the batch
   medians, the z-scaling and the decision threshold are all computed on the
   full dataset, including the pairs that are then scored. Work out which of
   these leak and by how much. Re-run with everything fit on a training split
   and applied to a held-out split, and report the drop.

3. **Claims 1 and 2 compare two different estimands across cohorts.**
   PXD001064 uses all cross-batch pairs; PXD004242 uses within-plate-group
   cross-plate pairs under one rotation assignment. Decide whether 0.9692 and
   0.9733 are legitimately comparable. If they are not, say what should be
   reported instead.

4. **Swaps are label permutations, not real swaps.** Does permuting labels
   produce the same difficulty as a physical mix-up? Think about what a real
   swap does that a permutation does not, and whether the simulation is
   optimistic, pessimistic, or neither.

5. **Donor identity is nested in plate group in PXD004242.** Verify that the
   within-group restriction actually removes the confound rather than hiding
   it, and that no residual plate information leaks into the score.

## ALSO CHECK

- Re-derive claims 1–6 from the raw inputs yourself. Do not accept the frozen
  JSON summaries in data/; they were written by the code you are auditing.
- `export_pxd004242_pairs.py` asserts the AUCs it should produce. Those frozen
  values have never been independently verified. Verify them.
- The MZ/DZ/unrelated ladder (claim 3): is the relatedness labelling correct,
  and are the DZ and unrelated classes properly disjoint?
- Figure 2 standardises each panel to its own null distribution. Is that
  defensible, or does it flatter the comparison? Would a reader be misled?
- Read `proteotrace/core.py` and `proteotrace/report.py` as code. Look for
  off-by-one errors, silent NaN handling, pairs counted twice, and anything
  where a verdict depends on row order.
- Run the tool on a cohort where you know the answer and the author did not:
  construct one, or degrade an existing one, and see whether it behaves.
- Check the figures against the numbers they claim to show. Open the PNGs.

## WHAT I ALREADY KNOW IS WRONG — do not spend time rediscovering, but do check that the stated limits are the real ones

- Monozygotic twins are not separable.
- Single-run donors cannot be assessed.
- Both cohorts are plasma, from two labs.
- A rejected Bonferroni-style rival adjustment halved sensitivity; it is off
  by default. Confirm that was the right call.

Do not take these as the complete list. If a stated limit is understated,
say so.

## RULES OF EVIDENCE

- Show the command and the output for every finding. A claim without a
  reproduction is a hypothesis, label it as one.
- Separate "I checked this and it holds" from "I did not check this". Say
  which parts of the repository you did not look at.
- If a number in the README does not reproduce, say so even if the gap is
  small. Especially then.
- If you think a claim is fine, say that plainly. An audit that finds nothing
  in a section is a result, not a failure.
- Do not propose new experiments as a substitute for auditing the existing
  ones. Scope creep is how audits avoid their job.

## OUTPUT

A single markdown report:

1. **Verdict per claim** — holds / holds with caveat / does not hold /
   not checked. One line each, with the evidence.
2. **Findings**, ordered by how much they change the paper. For each: what is
   wrong, how you demonstrated it, how much it moves the headline number, and
   the cheapest honest fix.
3. **What the paper must say that it currently does not.**
4. **What I could not check, and why.**

Be direct. If the central claim does not survive, lead with that.
