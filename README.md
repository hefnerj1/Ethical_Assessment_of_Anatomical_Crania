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

## Notes on data files

The `data/raw/` directory preserves the uploaded files unchanged. The `data/processed/` directory includes helper CSVs that use the filenames expected by the combined script when they could be created directly from the uploaded files:

- `MSUANAT.csv`
- `ANAT003.csv` derived from `MSU-Anatomical_data-cleaned.xlsx`
- `Anat_Imputed_Data.csv` copied from `Anat_Imputed_Data-WITHOUT_msu_ANATOMICAL_DATA.csv`
- `data.csv` derived from `MSU-Anatomical_data-cleaned.xlsx` for the FAMD workflow

The combined script also references refined reference files that were not included in the uploaded set:

- `Anat_Imputed_Refined_Data_Ind_Only.csv`
- `Anat_Imputed_Refined_Data.csv`
- `Anat_Imputed_Refined_Data_How_Only.csv`

Place those files in `data/processed/` or the project root before running the corresponding refined analyses.

## Suggested workflow

1. Open `anatomical-craniometric-analysis.Rproj` in RStudio.
2. Confirm required packages are installed.
3. Run `Anatomical-Craniometric-Combined.r` from the project root.
4. Save generated files in `outputs/` or allow the script to create root-level output files ignored by git.
