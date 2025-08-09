# 450k Kidney (KIRC) EWAS Template

This repository mirrors the F1000Research workflow by Maksimovic, Phipson & Oshlack (2016) and adapts it for **kidney renal clear cell carcinoma (KIRC)**. It can: 
- Clean and standardize your **IDAT** folder and **SampleSheet**,
- Run **QC → normalization → filtering**, and
- (Optionally) perform **DMP/DMR, GO, DV, and cell-type** steps.

## Quick start

1. Put your raw `.idat` files (any nested structure) inside: `data/idats_raw/`
2. Put or edit your phenotype in: `config/phenotype_kirc.csv`
   - This template already includes your 443 KIRC rows with columns like `Sentrix ID`, `Sentrix position`, `Type`.
3. In R (4.3+ recommended):
```r
source('run_all.R')    # CLEAN_ONLY = TRUE will run through filtering and save clean objects
```
4. Flip `CLEAN_ONLY <- FALSE` in `run_all.R` to run full analysis (DMPs, DMRs, GO, DV, cell-type notes).

## Outputs
- `outputs/qc/` – detection P plots, density, MDS, qcReport
- `outputs/normalized/` – `rgSet`, `mSet`, `beta`, `M` saved as RDS
- `outputs/results/` – tables and plots for DMPs/DMRs/GO/DV

## Notes for solid tumors
- **Cell-type composition** methods bundled with `minfi` are trained on blood; for solid tumors (kidney), reference-based estimates are not appropriate. We include placeholders and discussion.
- Sex chromosomes can be retained or filtered depending on your study question. We keep them by default and include a toggle.

