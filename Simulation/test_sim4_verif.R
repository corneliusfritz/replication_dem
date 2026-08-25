rm(list = ls())
source("Simulation/helper.R")
library(redeem)
library(peakRAM)
cluster_rds <- readRDS("Simulation/Results/simulation_4_2.RDS")
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

time_changepoints <- seq(0, 1000, length.out = 10)[-c(1, 10)]

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

cutpoints <- seq(from = 10000, to = 100000, length.out = 10)

res_one <- simulation_m(
  x = 1,
  cutpoints = cutpoints,
  formula_0_1 = formula_0_1,
  formula_1_0 = formula_1_0,
  coef_0_1 = coef_0_1,
  coef_1_0 = coef_1_0,
  popularity_0_1 = popularity_0_1,
  popularity_1_0 = popularity_1_0,
  n_nodes = n_nodes,
  nr = FALSE,
  simultaneous_interactions = TRUE,
  verbose = FALSE
)

diff_events <- max(abs(res_one$number_events - cluster_first_run$number_events))

cat(sprintf("Difference in event count arrays: %d\n", diff_events))

if (diff_events == 0) {
  cat("SUCCESS: Simulation 4 verification is successful and matches the cluster results perfectly!\n")
} else {
  stop("FAILURE: Simulated event cut sizes do not match.")
}
