# Replication Materials: Scalable Durational Event Models: Application to Physical and Digital Interactions

## Overview

This repository contains the R scripts and package redeem to replicate all results from the respective paper. 

---

## 1. Prerequisites

Note, that we generally assume that the working directory is the main folder of the replication folder not the specific folders for the application and simulation. 
Some of the scripts are configured for high-performance computing environments with parallelization.
Therefore, the scripts (simulation_1_final.R, simulation_2_final.R, etc.) would be taking a long time to run on a desktop, since they make use of a served with 50 cores. 

## 2. Session Info 
> sessionInfo()
R version 4.5.0 (2025-04-11)
Platform: aarch64-apple-darwin20
Running under: macOS 26.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Dublin
tzcode source: internal

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cowplot_1.2.0      dplyr_1.1.4        peakRAM_1.0.2      gridExtra_2.3      MASS_7.3-65       
 [6] stringi_1.8.7      readr_2.1.6        patchwork_1.3.2    Hmisc_5.2-5        latex2exp_0.9.8   
[11] scales_1.4.0       ggExtra_0.11.0     ggpubr_0.6.2       ggplot2_4.0.1      stringr_1.6.0     
[16] redeem_1.1.0       RcppProgress_0.4.2 data.table_1.18.0 

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1   farver_2.1.2       S7_0.2.1           fastmap_1.2.0      promises_1.5.0    
 [6] digest_0.6.39      rpart_4.1.24       mime_0.13          lifecycle_1.0.5    cluster_2.1.8.1   
[11] Cairo_1.7-0        survival_3.8-3     magrittr_2.0.4     compiler_4.5.0     rlang_1.1.7       
[16] tools_4.5.0        yaml_2.3.12        knitr_1.51         ggsignif_0.6.4     labeling_0.4.3    
[21] htmlwidgets_1.6.4  RColorBrewer_1.1-3 abind_1.4-8        miniUI_0.1.2       withr_3.0.2       
[26] foreign_0.8-90     purrr_1.2.1        nnet_7.3-20        xtable_1.8-4       colorspace_2.1-2  
[31] isoband_0.3.0      cli_3.6.5          rmarkdown_2.30     ragg_1.5.0         generics_0.1.4    
[36] otel_0.2.0         rstudioapi_0.17.1  tzdb_0.5.0         splines_4.5.0      base64enc_0.1-3   
[41] vctrs_0.6.5        Matrix_1.7-4       carData_3.0-5      car_3.1-3          hms_1.1.4         
[46] rstatix_0.7.3      Formula_1.2-5      htmlTable_2.4.3    systemfonts_1.3.1  tidyr_1.3.2       
[51] glue_1.8.0         gtable_0.3.6       later_1.4.5        tibble_3.3.1       pillar_1.11.1     
[56] htmltools_0.5.9    R6_2.6.1           textshaping_1.0.4  shiny_1.12.1       evaluate_1.0.5    
[61] lattice_0.22-7     backports_1.5.0    broom_1.0.11       httpuv_1.6.16      Rcpp_1.1.1        
[66] nlme_3.1-168       checkmate_2.3.3    mgcv_1.9-4         xfun_0.55          pkgconfig_2.0.3   
