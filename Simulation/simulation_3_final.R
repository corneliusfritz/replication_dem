# Simulation Study 3: Computational Scalability (Execution Time & RAM vs Node Cohort Sizes N)
# This script compares the efficiency of the proposed three-step Blockwise estimator
# against the standard legacy Newton-Raphson (NR) algorithm as network size N increases.

rm(list=ls())
source("Simulation/helper.R")
library(stringr)
library(readr)
library(stringi)
library(MASS)
library(redeem)
library(parallel)
library(peakRAM)
library(ggpubr)
library(dplyr)
library(cowplot)
library(patchwork)

# Set random seed to ensure exact reproducibility
set.seed(123456)

# Grid of node sizes N to test scalability
n_nodes_poss <- c(50, 75, 100, 125, 150, 175)

simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE
estimate_newton <- FALSE
K <- 50             # Parallel threads
verbose <- FALSE

# Setup fork cluster and load redeem package across threads
clust <- makeForkCluster(K, outfile = "outfile.txt")
clusterCall(clust, function() { library(redeem) })

# Loop through each network size grid point to benchmark performance
for(n_nodes in n_nodes_poss){
  cat("Starting ", n_nodes, "\n")
  
  # Measure time and RAM profiles under Blockwise and Newton-Raphson approaches
  run_simulation_for_n_node_time(n_nodes, n_nodes, simultaneous_interactions, it_max,
                                 accelerated, estimate_newton, K*20, verbose, clust, path = "Simulation/Results/")
  
  cat("Done with ", n_nodes, "\n")
}

# Cleanup cluster resources
file.remove("outfile.txt")
stopCluster(clust)

