# Runner script for Replication Revision with Results
# This master script sources the R packages installer, handles directory setup,
# manages computationally intensive simulation tasks, and executes the core 
# data processing, model estimation, out-of-sample validation, and plotting scripts.

# 0. Preparation
# Dynamically determine active script directory to resolve paths correctly
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


library(RhpcBLASctl)
RhpcBLASctl::blas_get_num_procs()
blas_set_num_threads(1)
omp_set_num_threads(1)
dirs_to_check <- c("Application/Plots", "Application/Results", "Simulation/Plots", "Simulation/Results")
for (d in dirs_to_check) {
  if (!dir.exists(d)) {
    message("Creating missing directory: ", d)
    dir.create(d, recursive = TRUE)
  }
}
source("install_packages.R")


# --- 1. Simulation folder ---
message("\n--- Starting Simulation scripts ---")

# Scenario 1: Basic DEM simulation across different node sizes (50 to 1000 nodes)
message("1/7 Sourcing: Simulation/simulation_1_final.R")
source("Simulation/simulation_1_final.R")

# Scenario 2: Simulated model performance comparison under varying parameters
message("2/7 Sourcing: Simulation/simulation_2_final.R")
source("Simulation/simulation_2_final.R")

# Scenario 3: Baseline comparisons and benchmarking of the Newton-Raphson vs Blockwise estimators
message("3/7 Sourcing: Simulation/simulation_3_final.R")
source("Simulation/simulation_3_final.R")

# Scenario 4: Alternative time-baseline simulation
message("4/7 Sourcing: Simulation/simulation_4_final.R")
source("Simulation/simulation_4_final.R")

# Scenario 5: Model selection criteria simulations (AIC / BIC accuracy)
message("5/7 Sourcing: Simulation/simulation_5_final.R")
source("Simulation/simulation_5_final.R")

# Simulation Plots: Plots the results of Scenario 2
message("6/7 Sourcing: Simulation/simulation_2_plot.R")
source("Simulation/simulation_2_plot.R")

# General Simulation Plots: Main simulation output visualizer
message("7/7 Sourcing: Simulation/simulation_plots.r")
source("Simulation/simulation_plots.r")


# --- 2. Application folder ---
message("\n--- Starting Application scripts ---")

message("1/4 Sourcing: Application/data_cleaning.R")
# Standardizes proximity Bluetooth scans, enforces undirected properties, and extracts continuous co-location events.
source("Application/data_cleaning.R")

message("2/4 Sourcing: Application/script_estimation.R")
# Formulates and fits the full Durational Event Models (DEM) for Calls and Proximity, and fits REM for SMS.
source("Application/script_estimation.R")

message("3/4 Sourcing: Application/script_oos.R")
# Performs out-of-sample (80% train, 20% test) splits and fits competing GoF baseline models.
source("Application/script_oos.R")

message("4/4 Sourcing: Application/script_plots.R")
# Generates all figures (degrees, trend lines, pop marginals) and outputs LaTeX results tables.
source("Application/script_plots.R")

message("\nAll scripts sourced successfully.")
message("Results saved in Application/Results and Simulation/Results.")
message("Plots saved in Application/Plots and Simulation/Plots.")

