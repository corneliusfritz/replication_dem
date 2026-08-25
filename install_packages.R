cran_packages <- c(
  "data.table",
  "stringr",
  "readr",
  "stringi",
  "dplyr",
  "Hmisc",
  "scales",
  "ggplot2",
  "ggpubr",
  "ggExtra",
  "cowplot",
  "patchwork",
  "gridExtra",
  "Cairo",
  "MASS",
  "latex2exp",
  "igraph",
  "peakRAM", 
  "redeem"
)

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    } else {
      message(paste("Package already installed:", pkg))
    }
  }
}


install_if_missing(cran_packages)

