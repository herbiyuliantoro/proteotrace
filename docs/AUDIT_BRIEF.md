# Independent audit scope

This document records the scope of the independent audit completed on
2026-10-02, before the preprint was prepared. The audit examined the analysis
state at commit `5f5d5b4` and the ProteoTrace 0.1.0 implementation. Later fixes
are described in [`ANALYSIS_LOG.md`](ANALYSIS_LOG.md); findings are in
[`AUDIT_REPORT.md`](AUDIT_REPORT.md).

## Questions examined

- Do the reported same-donor versus different-donor results reproduce from the
  cohort data?
- Do feature selection, preprocessing, or threshold estimation use information
  from the pairs being scored?
- Are the simulated swaps representative of the stated use, and are affected
  runs counted correctly?
- Do the analysis code and figures match the methods and reported results?

## Reproduction

The scripts in [`../audit/`](../audit/) reimplemented the primary scoring
analysis from the PXD001064 matrix and sample map, and from the PXD004242
MaxQuant export and SDRF. The independent scoring code does not import the
project's analysis scripts or read its frozen JSON summaries. Tests of tool
behaviour call the ProteoTrace package and independently account for the
simulated swaps.

From the repository root, run:

```sh
cd audit
sh run_all.sh
```

The saved command outputs are in [`../audit/outputs/`](../audit/outputs/).
The report covers the versions listed above. It does not assess every later
change to the repository.
