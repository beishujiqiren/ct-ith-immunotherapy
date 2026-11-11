# Install the R packages used by the supplied analysis scripts as needed.
cran_packages <- c(
  "tidyverse", "caret", "pheatmap", "ggpubr", "openxlsx", "vegan",
  "Rtsne", "ggstatsplot", "patchwork", "scales", "oro.nifti",
  "ggnewscale", "readxl", "broom", "stringr", "zoo", "tibble",
  "pROC", "forestploter", "rms", "randomForest", "readr",
  "gtsummary", "writexl", "imputeTS", "mice", "tidyr", "car",
  "flextable", "officer", "logistf", "Seurat", "clusterProfiler",
  "CytoTRACE2", "slingshot", "CellChat"
)
missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ConsensusClusterPlus")
}
