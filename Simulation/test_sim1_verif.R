rm(list = ls())
source("Simulation/helper.R")
library(redeem)
cluster_rds <- readRDS("Simulation/Results/simulation_1_50.RDS")
cluster_first_run <- cluster_rds[[1]][[1]] # Retrieve the first parallel iteration run
n_nodes <- 50
set.seed(50)

continuous_cov <- rnorm(n = n_nodes)
continuous_cov <- outer(X = continuous_cov, Y = continuous_cov, FUN = function(x, y) {
  abs(x - y)
})
categorical_cov <- sample(x = 1:3, size = n_nodes, replace = TRUE)
categorical_cov <- outer(X = categorical_cov, Y = categorical_cov, FUN = function(x, y) {
  x == y
})

popularity_0_1 <- rnorm(n = n_nodes, sd = 1) - (6 + 0.1 * log(n_nodes))
popularity_1_0 <- rnorm(n = n_nodes, sd = 1) - (1.6 + 0.1 * log(n_nodes))

time_changepoints <- seq(0, 10000, length.out = 10)[-c(1, 10)]

baseline_0_1_gt <- -seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
baseline_1_0_gt <- seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]

formula_0_1 <- ~ degree + current_common_partners(transformation = "log") +
  dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov) + baseline(changepoints = time_changepoints)
coef_0_1 <- c(
  "current_common_partners" = -0.5,
  "continuous_cov" = 1,
  "categorical_cov" = 0.5
)
formula_1_0 <- ~ degree + inertia(transformation = "log") +
  dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov) + baseline(changepoints = time_changepoints)
coef_1_0 <- c(
  "inertia" = 0.5,
  "continuous_cov" = 0.5, "categorical_cov" = 0.5
)

res_one <- simulation_run(
  x = 1,
  formula_0_1 = formula_0_1,
  formula_1_0 = formula_1_0,
  coef_0_1 = coef_0_1,
  coef_1_0 = coef_1_0,
  baseline_0_1_gt = baseline_0_1_gt,
  baseline_1_0_gt = baseline_1_0_gt,
  n_nodes = n_nodes,
  popularity_0_1 = popularity_0_1,
  popularity_1_0 = popularity_1_0,
  estimate = FALSE,
  simultaneous_interactions = TRUE,
  verbose = FALSE,
  it_max = 1500,
  accelerated = FALSE
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
  cat("SUCCESS: Simulation 1 verification is successful and matches the cluster results perfectly!\n")
} else {
  stop("FAILURE: Coefficients do not match.")
}

