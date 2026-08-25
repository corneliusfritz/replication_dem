# Replication Materials: Simulation Studies for Scalable Durational Event Models - Simulation

## 1. Prerequisites

### Software
* **R** (Version 4.0.0 or higher recommended)

### Required R Packages
Ensure the following packages are installed. The `redeem` package is required for the core model estimation.

```r
install.packages(c("data.table", "stringr", "readr", "stringi", "MASS",
                   "parallel", "peakRAM", "ggpubr", "dplyr", "cowplot",
                   "patchwork", "scales", "latex2exp", "Hmisc", "Cairo",
                   "grid", "gridExtra"))
# Ensure the 'redeem' package is installed
library(redeem)
```

## 2. Replication Steps

### Simulation 1: Consistency and Bias (Varying N)
Script: simulation_1_final.R
Objective: Evaluates the consistency of the estimator as the number of nodes ($N$) increases from 50 to 700.
Output: Saves results as simulation_1_[N].RDS in Simulation/Results/.

### Simulation 2: Asymptotic Normality
Script: simulation_2_final.R
Objective: Verifies the asymptotic normality of the estimators on a fixed network size ($N=500$).
Output: Saves simulation_2.RDS in Simulation/Results/.
Script: simulation_2_plot.R
Objective: Plots the estimation error distributions and coverage probabilities from Simulation 2 and generate tha basis for Table 1 in the manuscript (note that some parts of the table were removed by hand).

### Simulation 3: Scalability (Time & Memory vs. N)
Script: simulation_3_final.R
Objective: Compares the computational efficiency (Execution Time and Peak RAM) of the proposed Blockwise estimation against a standard Newton-Raphson approach as the number of actors ($N$) increases.
Output: Saves simulation_3_[N].RDS in Simulation/Results/ and generates accuracy plots simulation_3_accuracy.png and efficiency plots simulation_3.png.

### Simulation 4: Scalability (Time & Memory vs. M)
Script: simulation_4_final.R
Objective: Assesses the scalability of the estimator as the number of events ($M$) increases, keeping the number of nodes fixed ($N=500$).
Output: Saves results to Simulation/Results/simulation_4_2.RDS and generates the visualization simulation_4.png (Time vs. Events and Memory vs. Events) in Simulation/Plots/.

### Simulation 5: Approximation Accuracy (Grid Cells)
Script: simulation_5_final.R
Objective: Assesses how the number of grid cells (steps of the baseline hazard approximation) affects coverage probability.
Output: Saves simulation_5.RDS and generates simulation_5.png in Simulation/Plots/.

### Combined Visualization
Script: simulation_plots.R
Objective: Aggregates results from Simulations 2 and 3 into Figure 2.
Output: Generates simulation_final.png containing panels for Estimation Error, Time Comparison, and Memory Comparison.

## 3. Computation Settings
The scripts are configured for high-performance computing environments with parallelization. 
Reduce K in the scripts (simulation_1_final.R, simulation_2_final.R, etc.) if running on a standard laptop.
In this version of the replication package, we also provide verfication scripts for the simulation results. These scripts are located in the Simulation/Verification folder and can be run to verify the results of the simulations (you basically run one iteration locally and compare it with the saved results). 
The verification scripts are: test_sim1_verif.R, test_sim2_verif.R, test_sim3_verif.R, test_sim4_verif.R, and test_sim5_verif.R.
