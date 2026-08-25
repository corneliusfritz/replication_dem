# Simulation Study 1: Consistency and Bias under Varying Node Sizes (N)
# This script evaluates the asymptotic consistency of the Durational Event Model (DEM)
# estimator as the network size N grows from 50 to 1000 nodes.
# It sets up parallel cluster threads and triggers the corresponding helper estimation routine.

rm(list = ls())
source("Simulation/helper.R")
library(stringr)
library(readr)
library(stringi)
library(MASS)
library(redeem)
library(parallel)

# Set random seed to ensure exact reproducibility
set.seed(12345)
# n_nodes_poss <- c(1000)
n_nodes_poss <- c(50, 100, 300, 500, 700, 1000)

simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE
estimate_newton <- FALSE
K <- 50 # Number of simulation iterations per node configuration
verbose <- FALSE

# Setup parallel processing using fork-based multi-threading for computational efficiency
clust <- makeForkCluster(K, outfile = "outfile.txt")
clusterCall(clust, function() {
  library(redeem)
})
n_nodes <- 50

# Iterate through each node size and estimate the DEM consistency stats
for (n_nodes in n_nodes_poss) {
  cat("Starting ", n_nodes, "\n")
  run_simulation_for_n_node(n_nodes, n_nodes, simultaneous_interactions, it_max,
    accelerated, estimate_newton, K, verbose, clust,
    path = "Simulation/Results/"
  )
  cat("Done with ", n_nodes, "\n")
}

# Cleanup cluster resources and temporary log output
file.remove("outfile.txt")
stopCluster(clust)
