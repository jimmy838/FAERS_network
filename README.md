# FAERS Network Analysis Code

This repository contains the numerical code used to estimate symptom co-occurrence and Bayesian networks for FAERS reports in the manuscript "Co-occurrence Patterns of Symptom and Functional Burden Associated with Tofacitinib and Upa". It includes no network plotting code. Interactive network display is provided separately through the R Shiny application.

## Analyses

- `scripts/01_cooccurrence_network.R` estimates Ising, Phi, and PPMI co-occurrence networks.
- `scripts/02_bayesian_network.R` estimates a bootstrap-averaged Bayesian network with `death` constrained as a downstream node.

Both scripts use relative paths. Set `ANALYSIS_DRUG` to `tofacitinib` or `upadacitinib`; set `NETWORK_PROJECT_ROOT` only when running the project from another working directory.

## Input

Place one prepared RDS file for each drug in `data/input/`:

```text
network_input_tofacitinib.rds
network_input_upadacitinib.rds
```

The required object is a named list with `results_pt`, `reac`, `drug_pids`, and `procedure_pts`. The exact column requirements are documented in `data/input/README.md`. Raw FAERS, CVARD, and MedDRA files are not included in this repository.

## Run

From the repository root:

```r
install.packages(c("bnlearn", "data.table", "dplyr", "igraph", "IsingFit",
                   "openxlsx", "tidyr"))

Sys.setenv(ANALYSIS_DRUG = "tofacitinib")
source("scripts/01_cooccurrence_network.R")
source("scripts/02_bayesian_network.R")
```

The scripts write edge lists, adjacency matrices, graph RDS files, node and cluster statistics, network density/modularity/small-worldness measures, inter-network edge and cluster agreement, and Bayesian-network summaries to `results/network/`. The bootstrap seed and iteration count are defined in `config.R`.

The full output-file list and Shiny upload instructions are in `data/output/README.md`.

## Scope

This is a compact analysis package prepared for reproducibility. Community labels and node sizes are stored as graph attributes for downstream use. Data preparation, signal detection, time-to-onset analysis, and interactive visualization are maintained separately. Layout coordinates are intentionally left to the Shiny display layer.

## Compact demo

`demo/network_bayesian_demo.R` follows the same analysis sequence with simulated binary event data. It is included only to show the workflow and does not reproduce the study results.
