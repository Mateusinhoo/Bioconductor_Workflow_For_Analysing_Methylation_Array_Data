# methyCleanr — Universal 450k/EPIC methylation cleaner + EWAS starter

**Goal:** a single repository that **anyone** can use to (a) *clean* raw Illumina methylation arrays (IDATs → QC/normalize/filter) and (b) optionally run the classic analysis steps (DMP, DMR, GO, DV). It also supports starting from an **existing `rgSet`** rather than raw IDATs.

This template implements the F1000Research workflow (Maksimovic–Phipson–Oshlack) in a config-first way and ships with a **Conda environment** for easy setup.

---

## Why this structure? (Inputs/Outputs design)

### Inputs (two modes)
- **Mode A — Raw arrays (recommended):**  
  - `data/idats_raw/` → place all `.idat` files here (recursive allowed)  
  - `config/phenotype.csv` → must include at least `Slide` and `Array` (or `Sentrix ID`, `Sentrix position`). Optional: `Sample_Name`, `Sample_Group`, plus any custom metadata.
- **Mode B — Existing object:**  
  - `config/rgSet.rds` → a previously created `RGChannelSet` (e.g., from another pipeline)  
  - `config/phenotype.csv` → **must** include `Sample_Name` matching the columns of the object and a grouping column if you plan to run DMP/DMR/DV.

> Rationale: Most users have IDATs; some already have an `rgSet`. Supporting both maximizes reusability.

### Outputs (clean + analysis-ready)
- `outputs/rgSet.rds` — raw-level RGChannelSet saved
- `outputs/mSet.rds` — normalized `MethylSet`
- `outputs/beta_filtered.rds`, `outputs/M_filtered.rds` — **clean** matrices ready for downstream work
- `outputs/detectionP.rds`, `outputs/filter_summary.csv`, `outputs/sample_sheet_used.csv`
- QC/plots: `outputs/qcReport.pdf` (best-effort), `beta_density.png`, `mds.png`, `mean_detection_p_per_sample.png`
- If analysis enabled: `outputs/DMP_results.csv`, `outputs/DMR_ranges.rds`, `outputs/GO_results.csv`, `outputs/DV_results.csv`

> Rationale: RDS is compact and preserves attributes; CSV is used only for tabular results. Matrices are not exported as CSV by default because 450k×N can be huge.

### Extensibility
- Optional extra filters (cross-reactive, SNPs) via config
- Optional export to CSV/Parquet (off by default) via config
- Modular: advanced users can split `run_all.R` into per-step scripts; we keep a single entrypoint for simplicity.

---

## Quick start (Conda)

```bash
# 1) Create & activate env (Conda or Mamba)
mamba env create -f environment.yml
conda activate methycleanr

# 2) Put your inputs
#   - IDATs under data/idats_raw/
#   - phenotype CSV path in config/config.yaml

# 3) Run
Rscript -e "source('run_all.R')"
```

## Configuration

Edit `config/config.yaml`. Key fields:

```yaml
project_name: "example_project"

input:
  mode: "idat"                 # "idat" or "rgset"
  idat_dir: "data/idats_raw"
  phenotype_csv: "config/phenotype.csv"
  rgset_rds: "config/rgSet.rds"     # used only if mode = rgset

array_type: "auto"             # "auto", "450k", or "epic"

qc:
  detP_threshold: 0.01
  detP_fraction_pass: 0.90
  max_na_fraction: 0.05
  remove_sex: false            # drop chrX/Y probes if true

extra_filters:
  cross_reactive_list: ""      # optional path to probe list (one CpG per line)
  snp_probes_list: ""          # optional path to probe list (one CpG per line)

steps:
  load: true
  qc: true
  normalize: true
  explore: true
  filter: true
  idat_sanitize: false         # copy/rename non-standard IDATs to a clean folder (optional)
  dmp: false
  dmr: false
  go: false
  dv: false
  celltype_note: true

model:
  group_column: "Sample_Group"
  contrast: "Case- Control"

export:
  beta_csv: false              # WARNING: Huge files for many samples
  M_csv: false
  parquet: false               # enable if you add r-arrow/arrow
```

---

## Files & folders

```
methycleanr/
├── README.md
├── LICENSE
├── environment.yml                 # Conda env for all required packages
├── .github/workflows/ci.yml        # CI smoke test
├── .gitignore
├── config/
│   ├── config.yaml                 # knobs for steps/thresholds and inputs
│   ├── phenotype.csv               # example; replace with your own
│   └── rgSet.rds                   # optional (if starting from an object)
├── data/
│   ├── idats_raw/                  # <— drop .idat files here
│   └── annotation/                 # optional: probe lists
├── outputs/                        # generated
├── R/                              # room for modular scripts (optional)
├── scripts/
│   └── check_env.R                 # prints sessionInfo for debugging
└── run_all.R                       # single entrypoint
```

---

## Notes on cell-type composition

Houseman-style deconvolution methods are trained on blood references and **are not appropriate for solid tumors** like kidney. The pipeline writes a note if enabled and leaves a hook for users to plug in organ-specific or reference-free methods.

---

## License

MIT — do whatever you want, attribution appreciated.


### Selecting a subset of samples (from huge mixed IDAT folders)
Use `config/config.yaml -> select` to choose exactly which samples to load **before** reading IDATs:
- by **filter** expression on your phenotype columns (recommended),
- by **row indices** (e.g., 1500–2000), or
- by explicit **IDs** (e.g., a list of TCGA barcodes).
Only selected samples will be loaded from disk; others are ignored, even if thousands of IDATs exist.


## Required input: `config/samples.csv` (can be made in Excel)
You **must** provide a sample sheet at `config/samples.csv`. You can create it in Excel or Google Sheets and then **File → Save As → CSV**. Minimum required columns:
- `Slide` — the Sentrix (Slide) ID, e.g., `6264509100`
- `Array` — the Sentrix position, e.g., `R01C01`

Optional but recommended:
- `Sample_Name` — your unique sample identifier (e.g., TCGA barcode). If omitted, it will be auto-filled.
- `Sample_Group` — the group label used for plots/contrasts (e.g., `Primary solid Tumor`, `Solid Tissue Normal`).

> If your sheet uses headers `Sentrix ID` and `Sentrix position`, the pipeline will automatically map them to `Slide` and `Array`.

### Default selection mode: **indices 1500–2000**
By default, `config/config.yaml` is set to **load only rows 1500 through 2000 of `config/samples.csv`** before accessing the IDATs. Change these numbers as needed:
```yaml
select:
  enable: true
  by: "indices"
  indices:
    start: 1500
    end: 2000
```
This is useful when you have a huge mixed folder of IDATs and only want a specific subset.
