# methyCleanr — universal IDAT cleaning + EWAS starter

**methyCleanr** is a config-driven template that implements the F1000Research methylation array workflow (Maksimovic, Phipson & Oshlack) in a way anyone can run on their own 450k/EPIC IDATs.

- Works with **any** project: drop your `.idat` files in `data/idats_raw/` and point `config/config.yaml` to your phenotype CSV.
- Runs **step-by-step**: Loading → QC → Normalisation → Exploration → Filtering → (optional) DMP/DMR/GO/DV/Cell-type.
- Ships with Docker and renv options for reproducibility.

## Quick start
1. Put your raw `.idat` files (recursive ok) under: `data/idats_raw/`
2. Edit `config/config.yaml`:
   - `phenotype_csv`: path to your sample sheet or phenotype table
   - tweak thresholds and steps
3. In R:
```r
source("run_all.R")         # uses config/config.yaml
```
4. Results go to `outputs/`.

## Phenotype file
- Preferred columns: `Sample_Name, Sample_Group, Slide, Array`.
- Minimum required: `Slide, Array` (we’ll auto-fill others).
- If you have TCGA/other metadata, you can include extra columns; they are carried through to targets.

## Why this approach?
- **Universal**: no kidney-specific assumptions; works for blood or solid tissues (with caveats noted for cell-type deconvolution).
- **Configurable**: turn steps on/off; set thresholds without editing code.
- **Extensible**: add your tissue-specific deconvolution or custom filters.

## License
MIT
