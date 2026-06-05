# Anatomical Craniometric Analysis

This repository contains the combined R workflow and source data files for the anatomical craniometric analysis.

## Repository structure

```text
.
├── Anatomical-Craniometric-Combined.r      # convenience copy of the combined script
├── R/                                      # analysis scripts
├── data/
│   ├── raw/                                # original uploaded source files
│   └── processed/                          # script-compatible CSVs available from the uploaded files
├── docs/                                   # project documentation
└── outputs/                                # generated outputs; ignored by git except .gitkeep
```

## R packages

The combined script uses the following R packages:

```r
install.packages(c(
  "HDMD", "mnormt", "psych", "MASS", "readr", "tibble", "fansi",
  "mice", "ggplot2", "ellipse", "FactoMineR", "factoextra", "missMDA"
))
```

