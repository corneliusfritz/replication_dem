# Simulation Study 2: Asymptotic Normality under Fixed Cohort Size (N = 500)
# This script evaluates whether the DEM parameter estimators converge asymptotically
# to a multivariate normal distribution under a fixed network size.
# It simulates network event lists and runs parallel model estimations to extract coefficients.

rm(list=ls())
source("Simulation/helper.R")
library(stringr)
library(readr)
library(stringi)
library(MASS)
library(redeem)
library(parallel)

# Set random seed to ensure exact reproducibility
set.seed(1234)

# --- 1. Simulation Inputs & Covariates Setup ---
n_nodes <- 500

# Continuous dyadic covariate: absolute differences between normal random actor features
continuous_cov = rnorm(n = n_nodes)
continuous_cov = outer(X = continuous_cov, Y = continuous_cov, FUN = function(x, y) { abs(x - y) })

# Categorical dyadic covariate: homophily matching across 3 discrete actor categories
categorical_cov = sample(x = 1:3, size = n_nodes, replace = TRUE)
categorical_cov = outer(X = categorical_cov, Y = categorical_cov, FUN = function(x, y) { x == y })

# Popularity (nodal degree random effects) for both incidence and dissolution processes
popularity_0_1 = rnorm(n = n_nodes, sd = 1) - (6 + 0.1 * log(n_nodes))
popularity_1_0 = rnorm(n = n_nodes, sd = 1) - (1.6 + 0.1 * log(n_nodes))

# --- 2. Ground-Truth Baseline Trends ---
# Baseline hazards changing at fixed time changepoints (intervals 1 to 10)
time_changepoints <- seq(0, 10000, length.out = 10)[-c(1, 10)]
baseline_0_1_gt <- -seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
baseline_1_0_gt <- seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]

# --- 3. Specify Ground-Truth Models & Coefficients ---
# formula_0_1: Incidence (starting an event) model
formula_0_1 =  ~ degree + current_common_partners(transformation = "log") +
  dyadic_cov(data = continuous_cov)  + dyadic_cov(data = categorical_cov)
coef_0_1= c("current_common_partners" = -0.5,
            "continuous_cov" = 1,
            "categorical_cov" = 0.5)

# formula_1_0: Dissolution (ending an active event) model
formula_1_0 =  ~ degree + number_interaction(transformation = "log") +
  dyadic_cov(data = continuous_cov)  + dyadic_cov(data = categorical_cov)
coef_1_0 = c("number_interaction" = 0.5,
             "continuous_cov" = 0.5, "categorical_cov" = 0.5)

# Save setting parameters for plot referencing
save.image(file = "Simulation/Results/setting_simulation_2.RData")

simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE
estimate_newton <- FALSE

K <- 100             # Parallel iterations
verbose <- FALSE

# --- 4. Parallel Cluster Execution ---
# Spin up fork cluster for parallel model estimation
clust <- makeForkCluster(K, outfile = "outfile.txt")
clusterCall(clust, function() { library(redeem) })

# Run simulation_model_selection on each parallel thread
res_simulation = parLapply(cl = clust, X = 1:(K*10), fun = simulation_model_selection,
                           formula_0_1 = formula_0_1, formula_1_0 = formula_1_0,
                           coef_0_1 = coef_0_1, coef_1_0 = coef_1_0,
                           time_changepoints = time_changepoints,
                           baseline_0_1_gt = baseline_0_1_gt,
                           baseline_1_0_gt = baseline_1_0_gt,
                           n_nodes = n_nodes, popularity_0_1 = popularity_0_1,
                           popularity_1_0 = popularity_1_0,
                           simultaneous_interactions = simultaneous_interactions,
                           verbose = FALSE, it_max = it_max, accelerated = accelerated)

# Save completed simulation RDS data
saveRDS(res_simulation, file = "Simulation/Results/simulation_2.RDS")

# Cleanup cluster resources
file.remove("outfile.txt")
stopCluster(clust)

