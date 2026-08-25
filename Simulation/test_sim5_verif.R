rm(list = ls())
source("Simulation/helper.R")
library(redeem)

# Load cluster RDS file
cluster_rds <- readRDS("Simulation/Results/simulation_5.RDS")
cluster_first_run <- cluster_rds[[1]] # Retrieve the first parallel iteration run
n_nodes <- 500
set.seed(1234)

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

formula_0_1 <- ~ degrees + current_common_partners(transformation = "log") +
  dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov)
coef_0_1 <- c(
  "current_common_partners" = -0.5,
  "continuous_cov" = 1,
  "categorical_cov" = 0.5
)
formula_1_0 <- ~ degrees + inertia(transformation = "log") +
  dyadic_cov(data = continuous_cov) + dyadic_cov(data = categorical_cov)
coef_1_0 <- c(
  "inertia" = 0.5,
  "continuous_cov" = 0.5, "categorical_cov" = 0.5
)

grid_cells <- round(exp(seq(from = 1, to = 7, length.out = 20)) + 1)

res_one <- simulation_baseline(
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
  simultaneous_interactions = TRUE,
  verbose = FALSE,
  it_max = 1000,
  accelerated = FALSE,
  grid_cells = grid_cells[1:3]
)

# Compare coefficients (Incidence and Duration) across the 3 grid cells
diff_0_1 <- max(abs(unlist(res_one$est_0_1_list) - unlist(cluster_first_run$est_0_1_list[1:3])))
diff_1_0 <- max(abs(unlist(res_one$est_1_0_list) - unlist(cluster_first_run$est_1_0_list[1:3])))

cat(sprintf("Difference in finite Incidence coefficients across baseline grids: %e\n", diff_0_1))
cat(sprintf("Difference in finite Duration coefficients across baseline grids: %e\n", diff_1_0))

if (diff_0_1 < 1e-7 && diff_1_0 < 1e-7) {
  cat("SUCCESS: Simulation 5 verification is successful and matches the cluster results perfectly!\n")
} else {
  stop("FAILURE: Coefficients do not match.")
}
