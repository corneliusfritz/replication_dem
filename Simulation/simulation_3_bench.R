
rm(list=ls())

source("Simulation/helper.R")
library(stringr)
library(readr)
library(stringi)
library(MASS)
library(redeem)
library(parallel)
library(bench)
library(ggpubr)
library(dplyr)
library(cowplot)
library(patchwork)
library(data.table)
library(Cairo)

set.seed(123456)


n_nodes_poss <- c(50, 75, 100, 125, 150, 175)

simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE
estimate_newton <- FALSE
K <- 50
verbose <- FALSE
clust <- makeForkCluster(K,outfile = "outfile.txt")
clusterCall(clust, function() {library(redeem); library(bench)})

for(n_nodes in n_nodes_poss){
  cat("Starting ",n_nodes,"\n")
  run_simulation_for_n_node_time_bench(n_nodes, n_nodes, simultaneous_interactions, it_max,
                            accelerated, estimate_newton, K*20, verbose, clust, path= "Simulation/Results/")
  cat("Done with ",n_nodes,"\n")
}

file.remove("outfile.txt")
stopCluster(clust)
