# Replication Materials: Scalable Durational Event Models - Application

## 1. Prerequisites

### Software
* **R** (Version 4.5 or higher recommended)

### Required R Packages
Ensure the following packages are installed before running the scripts. The analysis relies heavily on the `redeem` package for model estimation.

```r
install.packages(c("data.table", "stringr", "ggplot2", "ggpubr", "ggExtra", 
                   "scales", "latex2exp","redeem", "Hmisc", "patchwork", "Cairo", "igraph"))

# Ensure the 'redeem' package is installed
library(redeem)
```
## 2. Replication Steps

Execute the scripts in the following order to ensure dependencies (processed data and saved model objects) are available for subsequent steps.

### Step 1: Data Pre-processing

Script: data_cleaning.R
Function: Processes the raw Bluetooth symmetric data (bt_symmetric.csv).
Output: Generates Application/events.csv (processed proximity events) and Application/actor_data.csv (actor ID mapping).

### Step 2: Model Estimation
Script: script_estimation.R
Function: Fits the Durational Event Models to the Call, Proximity, and SMS datasets.
Output: Saves fitted model objects as .rds files 

### Step 3: Out-of-Sample Validation & Goodness of Fit
Script: script_oos.R
Function: Performs goodness-of-fit assessments and out-of-sample prediction tasks.
Details:Splits events into training (80%) and test sets.Fits comparative models for validation.
Output: Saves Goodness-of-Fit results as .rds files in Application/Results/ (e.g., gof_results_denmark.rds).

### Step 4: Visualization and Tables
Script: script_plots.R
Function: Generates the final figures and LaTeX tables used in the paper.
Output: All plots 

Depending on your hardware, particularly Step 2 and 3 may take 2-3 hours to run.

