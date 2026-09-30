# Compact demonstration of the network workflow.

source(file.path("R", "network_functions.R"))
library(bnlearn)
library(IsingFit)

set.seed(42)
n_reports <- 300L
syndrome <- rbinom(n_reports, 1, 0.35)
make_event <- function(p1, p0) rbinom(n_reports, 1, ifelse(syndrome == 1, p1, p0))
events <- cbind(
  pain = make_event(0.70, 0.15),
  fatigue = make_event(0.60, 0.12),
  falls = make_event(0.45, 0.08),
  limitation = make_event(0.50, 0.10),
  headache = rbinom(n_reports, 1, 0.18),
  nausea = rbinom(n_reports, 1, 0.16),
  infection = rbinom(n_reports, 1, 0.12),
  death = rbinom(n_reports, 1, 0.04)
)

ising <- IsingFit(events, progressbar = FALSE, plot = FALSE)$weiadj
phi <- calculate_phi_network(events)
ppmi <- calculate_ppmi_network_bootstrap(events, n_bootstrap = 100)

blacklist <- data.frame(
  from = "death",
  to = setdiff(colnames(events), "death")
)
bn_data <- as.data.frame(events)
bn_data[] <- lapply(bn_data, factor, levels = c(0, 1))
bn_strength <- boot.strength(
  bn_data, R = 100, algorithm = "hc",
  algorithm.args = list(score = "bic", blacklist = blacklist)
)
bn_average <- averaged.network(bn_strength)

cat("Ising edges:", sum(ising != 0) / 2, "\n")
cat("Phi edges:", sum(phi != 0) / 2, "\n")
cat("PPMI edges:", sum(ppmi != 0) / 2, "\n")
cat("Bayesian edges:", nrow(arcs(bn_average)), "\n")
