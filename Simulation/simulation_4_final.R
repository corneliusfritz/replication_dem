# Simulation Study 4: Computational Scalability (Execution Time & RAM vs Event Counts M)
# This script assesses the scalability of the proposed Blockwise estimator as the number 
# of event records (M) increases from 10k to 100k, under a fixed node cohort size (N = 500).
# It profiles CPU time and Peak RAM using peakRAM and compiles boxplots of the results.

rm(list=ls())
source("Simulation/helper.R")
library(stringr)
library(readr)
library(stringi)
library(MASS)
library(redeem)
library(parallel)
library(ragg)
library(ggpubr)
library(scales)
library(patchwork)
library(data.table)

# Set random seed to ensure exact reproducibility
set.seed(12345)
seed <- 1234
simultaneous_interactions <- TRUE
it_max <- 1500
accelerated <- FALSE
estimate_newton <- FALSE
verbose <- FALSE
n_nodes <- 500

cat("Starting ", n_nodes, "\n")
set.seed(seed)

# --- 1. Set up Covariates ---
continuous_cov = rnorm(n = n_nodes)
continuous_cov = outer(X = continuous_cov, Y = continuous_cov, FUN = function(x, y) { abs(x - y) })
categorical_cov = sample(x = 1:3, size = n_nodes, replace = TRUE)
categorical_cov = outer(X = categorical_cov, Y = categorical_cov, FUN = function(x, y) { x == y })

popularity_0_1 = rnorm(n = n_nodes, sd = 1) - (6 + 0.1 * log(n_nodes))
popularity_1_0 = rnorm(n = n_nodes, sd = 1) - (1.6 + 0.1 * log(n_nodes))

time_changepoints <- seq(0, 1000, length.out = 10)[-c(1, 10)]

baseline_0_1_gt <- -seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]
baseline_1_0_gt <- seq(from = 0, to = -0.1, length.out = length(time_changepoints) + 1)[-1]

# Specify ground truth formula terms
formula_0_1 =  ~ degrees + current_common_partners(transformation = "log") +
  dyadic_cov(data = continuous_cov)  + dyadic_cov(data = categorical_cov)
coef_0_1= c("current_common_partners" = -0.5,
            "continuous_cov" = 1,
            "categorical_cov" = 0.5)

formula_1_0 =  ~ degrees +  inertia(transformation = "log") +
  dyadic_cov(data = continuous_cov)  + dyadic_cov(data = categorical_cov)
coef_1_0 = c("inertia" = 0.5,
             "continuous_cov" = 0.5, "categorical_cov" = 0.5)

# --- 2. Define the Event Grid (10k to 100k) ---
cutpoints <- seq(from = 10000, to = 100000, length.out = 10)
verbose <- TRUE
K <- 100

# --- 3. Parallel Cluster Execution ---
clust <- makeForkCluster(K, outfile = "outfile.txt")
clusterCall(clust, function() { library(redeem) })

# Run simulation_m to profile Blockwise estimation across increasing event cuts
res <- parLapply(cl = clust, X = 1:(K*20), fun = simulation_m, cutpoints = cutpoints,
                 formula_0_1 = formula_0_1, formula_1_0 = formula_1_0,
                 coef_0_1 = coef_0_1, coef_1_0 = coef_1_0, nr = FALSE,
                 n_nodes = n_nodes, popularity_0_1 = popularity_0_1,
                 popularity_1_0 = popularity_1_0, verbose = FALSE,
                 simultaneous_interactions = simultaneous_interactions)

# Save scalability RDS outputs
saveRDS(res, file = "Simulation/Results/simulation_4_2.RDS")
stopCluster(clust)

# --- 4. Plotting & Compilation of Boxplots ---
res <- readRDS(file = "Simulation/Results/simulation_4_2.RDS")
res <- mapply(x = res, y = seq_along(res), function(x, y) { x$id = y; return(x) }, SIMPLIFY = FALSE)
res <- do.call("rbind", res)
vals <- seq(10000, 100000, by = 10000)
labels_k <- paste0(vals / 1000, "k")
res <- as.data.table(res)

# Boxplot of CPU time needed vs Event Count
a <- ggplot(res, aes(y = time_needed, x = factor(number_events))) +
  geom_boxplot() +
  theme_pubr(base_size = 17) +
  xlab("Number of events (M)") +
  ylab("Time in seconds") +
  scale_x_discrete(breaks = seq(10000, 100000, by = 10000), label = labels_k)

# Boxplot of Peak Memory used vs Event Count
b <- ggplot(res, aes(y = peak_mem, x = factor(number_events))) +
  geom_boxplot() +
  theme_pubr(base_size = 17) +
  xlab("Number of events (M)") +
  ylab("Peak memory used in MB") +
  scale_x_discrete(breaks = seq(10000, 100000, by = 10000), label = labels_k)

# Combine both boxplots side-by-side and save to file
ggsave("Simulation/Plots/simulation_4.png", a | b, width = 14, height = 5)

