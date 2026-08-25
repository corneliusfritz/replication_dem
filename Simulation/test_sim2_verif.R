
rm(list = ls())
source("Simulation/helper.R")
library(redeem)

cluster_rds <- readRDS("Simulation/Results/simulation_2.RDS")
cluster_first_run <- cluster_rds[[1]] # Retrieve the first parallel iteration run

load(file = "Simulation/Results/setting_simulation_2.RData")
simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE

# Run exactly 1 iteration (x = 1) locally
res_one <- simulation_model_selection(
  x = 1,
  formula_0_1 = formula_0_1,
  formula_1_0 = formula_1_0,
  coef_0_1 = coef_0_1,
  coef_1_0 = coef_1_0,
  time_changepoints = time_changepoints,
  baseline_0_1_gt = baseline_0_1_gt,
  baseline_1_0_gt = baseline_1_0_gt,
  n_nodes = n_nodes,
  popularity_0_1 = popularity_0_1,
  popularity_1_0 = popularity_1_0,
  simultaneous_interactions = simultaneous_interactions,
  verbose = FALSE,
  it_max = it_max,
  accelerated = accelerated
)

finite_idx_0_1 <- which(is.finite(res_one$coef_0_1_core))
finite_idx_1_0 <- which(is.finite(res_one$coef_1_0_core))

inf_match_0_1 <- identical(is.infinite(res_one$coef_0_1_core), is.infinite(cluster_first_run$coef_0_1_core))
inf_match_1_0 <- identical(is.infinite(res_one$coef_1_0_core), is.infinite(cluster_first_run$coef_1_0_core))

diff_0_1 <- max(abs(res_one$coef_0_1_core[finite_idx_0_1] - cluster_first_run$coef_0_1_core[finite_idx_0_1]))
diff_1_0 <- max(abs(res_one$coef_1_0_core[finite_idx_1_0] - cluster_first_run$coef_1_0_core[finite_idx_1_0]))

cat(sprintf("Difference in finite Incidence coefficients: %e\n", diff_0_1))
cat(sprintf("Difference in finite Duration coefficients: %e\n", diff_1_0))
cat(sprintf("Identical Infinite coefficient structures: %s\n", inf_match_0_1 && inf_match_1_0))

if (inf_match_0_1 && inf_match_1_0 && diff_0_1 < 1e-7 && diff_1_0 < 1e-7) {
  cat("SUCCESS: Simulation 2 verification is successful and matches the cluster results perfectly!\n")
} else {
  stop("FAILURE: Coefficients do not match.")
}

