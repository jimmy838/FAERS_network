# Network input contract

Each file must be an RDS containing a named list with these elements:

```r
list(
  results_pt = data.frame(
    event = character(),
    IC_lower = numeric(),
    ROR_lower = numeric(),
    PRR_median = numeric(),
    chi_square = numeric(),
    EB05 = numeric()
  ),
  reac = data.frame(
    primaryid = character(),
    pt = character()
  ),
  drug_pids = character(),
  procedure_pts = character()
)
```

`results_pt` contains one row per preferred term and the signal-detection fields used by the network filter. `reac` contains one row per report and preferred term. `drug_pids` identifies reports for the study drug, and `procedure_pts` lists preferred terms removed before network estimation.

The files are intentionally excluded from version control because they may contain report-level pharmacovigilance data.
