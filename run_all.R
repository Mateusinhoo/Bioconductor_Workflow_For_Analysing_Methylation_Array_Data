# run_all.R — methyCleanr (universal; with sample selection)
suppressPackageStartupMessages({
  library(yaml); library(minfi); library(limma); library(DMRcate); library(missMethyl)
  library(ggplot2); library(matrixStats); library(R.utils); library(optparse)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# -------- CLI (optional) --------
opt_list <- list(
  make_option(c("-c","--config"), type="character", default="config/config.yaml",
              help="Path to config.yaml [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
cfg_path <- opt$config

# -------- Load config --------
cfg <- yaml::read_yaml(cfg_path)

prj        <- cfg$project_name %||% "project"
input      <- cfg$input
steps      <- cfg$steps
qc         <- cfg$qc
extra      <- cfg$extra_filters
export_opt <- cfg$export
array_type <- tolower(cfg$array_type %||% "auto")
select     <- cfg$select

dir.create("outputs", showWarnings = FALSE)

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n")

# -------- Helpers --------
normalize_pheno <- function(df) {
  n <- names(df)
  n <- sub("^Sentrix ID$","Slide", n, ignore.case = TRUE)
  n <- sub("^Sentrix position$","Array", n, ignore.case = TRUE)
  names(df) <- n
  df
}
require_cols <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop("Missing columns in phenotype CSV: ", paste(miss, collapse=", "))
}
apply_selection <- function(pheno, sel) {
  if (is.null(sel) || !isTRUE(sel$enable)) return(pheno)
  by <- tolower(sel$by %||% "filter")
  if (by == "filter") {
    filt <- sel$filter %||% ""
    if (nzchar(filt)) {
      keep <- with(pheno, eval(parse(text = filt)))
      keep[is.na(keep)] <- FALSE
      pheno <- pheno[keep, , drop = FALSE]
    }
  } else if (by == "indices") {
    st <- as.integer(sel$indices$start %||% 1)
    en <- as.integer(sel$indices$end %||% nrow(pheno))
    st <- max(1, st); en <- min(nrow(pheno), en)
    if (st <= en) pheno <- pheno[st:en, , drop = FALSE] else pheno <- pheno[0,]
  } else if (by == "ids") {
    idcol <- sel$ids_column %||% stop("ids_column must be set for select$by='ids'")
    ids <- character()
    if (!is.null(sel$ids_file) && nzchar(sel$ids_file) && file.exists(sel$ids_file)) {
      ids <- unique(trimws(readLines(sel$ids_file)))
    }
    if (!is.null(sel$ids_list) && length(sel$ids_list) > 0) {
      ids <- unique(c(ids, unlist(sel$ids_list)))
    }
    if (length(ids) == 0) stop("No IDs provided for selection.")
    if (!idcol %in% names(pheno)) stop("ids_column '", idcol, "' not found in phenotype.")
    pheno <- pheno[pheno[[idcol]] %in% ids, , drop = FALSE]
  } else {
    stop("Unknown selection method: ", by)
  }
  # Optional: keep certain types
  if (!is.null(sel$keep_types) && length(sel$keep_types) > 0 && "Type" %in% names(pheno)) {
    pheno <- pheno[pheno$Type %in% sel$keep_types, , drop = FALSE]
  }
  # Optional: limit N
  if (!is.null(sel$limit_n)) {
    n <- as.integer(sel$limit_n)
    if (!is.na(n) && n > 0 && n < nrow(pheno)) {
      pheno <- pheno[seq_len(n), , drop = FALSE]
    }
  }
  pheno
}

# -------- Inputs --------
mode <- tolower(input$mode %||% "idat")
pheno <- read.csv(input$phenotype_csv %||% "config/phenotype.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
pheno <- normalize_pheno(pheno)

# Apply selection BEFORE building targets
pheno <- apply_selection(pheno, select)

if (mode == "idat") {
  require_cols(pheno, c("Slide","Array"))
  if (!"Sample_Name" %in% names(pheno)) pheno$Sample_Name <- seq_len(nrow(pheno))
  if (!"Sample_Group" %in% names(pheno)) pheno$Sample_Group <- "Group"
  targets <- pheno[, intersect(c("Sample_Name","Sample_Group","Slide","Array"), names(pheno)), drop = FALSE]
  write.csv(targets, "outputs/sample_sheet_used.csv", row.names = FALSE)
  idat_dir <- input$idat_dir %||% "data/idats_raw"
  
  if (isTRUE(steps$idat_sanitize)) {
    out <- "data/idats_clean"
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    files <- list.files(idat_dir, pattern = "(?i)idat(\\.gz)?$", full.names = TRUE, recursive = TRUE)
    for (f in files) {
      bn <- basename(f)
      bn <- sub("\\.IDAT(\\.GZ)?$", ".idat\\1", bn, ignore.case = TRUE)
      if (grepl("_R\\.", bn, ignore.case=TRUE)) bn <- sub("_R\\.", "_Red.", bn)
      if (grepl("_G\\.", bn, ignore.case=TRUE)) bn <- sub("_G\\.", "_Grn.", bn)
      file.copy(f, file.path(out, bn), overwrite = FALSE)
    }
    idat_dir <- out
  }
  
  logmsg("Reading IDATs from %s ...", idat_dir)
  if (!dir.exists(idat_dir)) stop("IDAT directory does not exist: ", idat_dir)
  rgSet <- read.metharray.exp(base = idat_dir, targets = targets, extended = TRUE, force = TRUE)
  saveRDS(rgSet, file = "outputs/rgSet.rds")
  
} else if (mode == "rgset") {
  rg_path <- input$rgset_rds %||% "config/rgSet.rds"
  if (!file.exists(rg_path)) stop("rgSet.rds not found at: ", rg_path)
  rgSet <- readRDS(rg_path)
  require_cols(pheno, c("Sample_Name"))
  if (!"Sample_Group" %in% names(pheno)) pheno$Sample_Group <- "Group"
  sm <- sampleNames(rgSet)
  rownames(pheno) <- pheno$Sample_Name
  if (!all(sm %in% rownames(pheno))) {
    stop("Phenotype CSV must contain all sample names present in rgSet (Sample_Name column).")
  }
  pheno <- pheno[sm, , drop=FALSE]
  targets <- pheno[, intersect(c("Sample_Name","Sample_Group","Slide","Array"), names(pheno)), drop = FALSE]
  write.csv(targets, "outputs/sample_sheet_used.csv", row.names = FALSE)
} else {
  stop("Unknown input mode: ", mode)
}

# -------- Array type detection --------
ann_choice <- NULL
if (tolower(array_type) == "auto") {
  ann <- annotation(rgSet)
  ann_choice <- if (grepl("EPIC", paste(ann, collapse=" "), ignore.case = TRUE)) "EPIC" else "450K"
} else {
  ann_choice <- toupper(array_type)
}
logmsg("Array type selected: %s", ann_choice)

ann_obj <- switch(ann_choice,
  "EPIC" = IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
  "450K" = IlluminaHumanMethylation450kanno.ilmn12.hg19,
  stop("Unknown array type: ", ann_choice)
)

# -------- QC --------
if (isTRUE(steps$qc)) {
  logmsg("Computing detection P-values...")
  detP <- detectionP(rgSet)
  saveRDS(detP, "outputs/detectionP.rds")
  png("outputs/mean_detection_p_per_sample.png", width = 1200, height = 700)
  barplot(colMeans(detP, na.rm = TRUE), las = 2, ylab = "Mean detection P-value")
  abline(h = qc$detP_threshold %||% 0.01, lty = 2)
  dev.off()
  try(qcReport(rgSet, sampNames = targets$Sample_Name, sampGroups = targets$Sample_Group,
               pdf = "outputs/qcReport.pdf"), silent = TRUE)
}

# -------- Normalization --------
if (isTRUE(steps$normalize)) {
  logmsg("Normalizing (Noob -> Quantile)...")
  mSet.noob <- preprocessNoob(rgSet)
  mSet <- preprocessQuantile(mSet.noob)
  saveRDS(mSet, "outputs/mSet.rds")
  beta <- getBeta(mSet); M <- getM(mSet)
  saveRDS(beta, "outputs/beta_raw.rds"); saveRDS(M, "outputs/M_raw.rds")
} else {
  mSet <- readRDS("outputs/mSet.rds")
  beta <- readRDS("outputs/beta_raw.rds"); M <- readRDS("outputs/M_raw.rds")
}

# -------- Explore --------
if (isTRUE(steps$explore)) {
  logmsg("Generating density & MDS plots...")
  png("outputs/beta_density.png", width = 1200, height = 700)
  plotDensities(beta, main = "Beta densities (post-normalization)")
  dev.off()
  png("outputs/mds.png", width = 1200, height = 700)
  plotMDS(getM(mSet), top = 1000, gene.selection = "common",
          labels = targets$Sample_Name, col = as.numeric(factor(targets$Sample_Group)))
  dev.off()
}

# -------- Filtering --------
if (isTRUE(steps$filter)) {
  logmsg("Filtering probes...")
  detP <- if (exists("detP")) detP else readRDS("outputs/detectionP.rds")
  keep_det <- rowMeans(detP < (qc$detP_threshold %||% 0.01), na.rm = TRUE) >= (qc$detP_fraction_pass %||% 0.90)
  keep_nonNA <- rowMeans(is.na(beta)) <= (qc$max_na_fraction %||% 0.05)
  ann <- getAnnotation(ann_obj)
  sex_mask <- ann$chr %in% c("chrX","chrY")
  keep_sex <- if (isTRUE(qc$remove_sex)) !sex_mask else rep(TRUE, length(sex_mask))
  idx <- match(rownames(beta), rownames(ann))
  keep <- keep_det & keep_nonNA & keep_sex[idx]
  
  if ((extra$cross_reactive_list %||% "") != "") {
    xr <- readLines(extra$cross_reactive_list)
    keep <- keep & !(rownames(beta) %in% xr)
  }
  if ((extra$snp_probes_list %||% "") != "") {
    sp <- readLines(extra$snp_probes_list)
    keep <- keep & !(rownames(beta) %in% sp)
  }
  
  beta_f <- beta[keep,,drop=FALSE]; M_f <- M[keep,,drop=FALSE]
  saveRDS(beta_f, "outputs/beta_filtered.rds")
  saveRDS(M_f, "outputs/M_filtered.rds")
  write.csv(data.frame(probes_kept=sum(keep)), "outputs/filter_summary.csv", row.names = FALSE)
  
  if (isTRUE(export_opt$beta_csv)) write.csv(beta_f, gzfile("outputs/beta_filtered.csv.gz"))
  if (isTRUE(export_opt$M_csv)) write.csv(M_f, gzfile("outputs/M_filtered.csv.gz"))
}

# -------- Modeling helpers --------
get_design <- function(targets, group_col) {
  g <- factor(targets[[group_col]])
  design <- model.matrix(~ 0 + g); colnames(design) <- levels(g)
  list(g=g, design=design)
}
parse_contrast <- function(contrast, design) {
  makeContrasts(contrasts = contrast, levels = design)
}

# -------- DMP --------
if (isTRUE(steps$dmp)) {
  logmsg("Probe-wise DMP (limma)...")
  Mwork <- if (isTRUE(steps$filter)) readRDS("outputs/M_filtered.rds") else readRDS("outputs/M_raw.rds")
  dd <- get_design(targets, cfg$model$group_column %||% "Sample_Group")
  cont <- parse_contrast(cfg$model$contrast %||% "Case- Control", dd$design)
  fit <- lmFit(Mwork, dd$design)
  fit2 <- contrasts.fit(fit, cont); fit2 <- eBayes(fit2)
  dmp <- topTable(fit2, number = Inf, adjust.method = "BH")
  write.csv(dmp, "outputs/DMP_results.csv")
  png("outputs/volcano_dmp.png", width = 1100, height = 800)
  with(dmp, { plot(logFC, -log10(P.Value), pch=20, xlab="logFC (M)", ylab="-log10 P"); abline(v=c(-1,1), lty=2); abline(h=1.3, lty=2) })
  dev.off()
}

# -------- DMR --------
if (isTRUE(steps$dmr)) {
  logmsg("Regional DMRs (DMRcate)...")
  Mwork <- if (isTRUE(steps$filter)) readRDS("outputs/M_filtered.rds") else readRDS("outputs/M_raw.rds")
  dd <- get_design(targets, cfg$model$group_column %||% "Sample_Group")
  cont <- parse_contrast(cfg$model$contrast %||% "Case- Control", dd$design)
  ann <- getAnnotation(ann_obj)
  myAnnot <- cpg.annotate(object="array", datatype="array", what="M", analysis.type="differential",
                          x=Mwork, design=dd$design, contrasts=TRUE, cont.matrix=cont, coef=1,
                          fdr=0.05, annotation=ann)
  dmrcoutput <- dmrcate(myAnnot, lambda=1000, C=2)
  ranges <- extractRanges(dmrcoutput, genome="hg19")
  saveRDS(ranges, "outputs/DMR_ranges.rds")
}

# -------- GO --------
if (isTRUE(steps$go)) {
  logmsg("GO testing (gometh)...")
  dmp <- read.csv("outputs/DMP_results.csv", stringsAsFactors = FALSE)
  Mwork <- if (isTRUE(steps$filter)) readRDS("outputs/M_filtered.rds") else readRDS("outputs/M_raw.rds")
  sig <- dmp$X[dmp$adj.P.Val < 0.05]
  go <- gometh(sig.cpg = sig, all.cpg = rownames(Mwork), collection="GO",
               array.type = ifelse(ann_choice=='EPIC','EPIC','450K'), plot.bias=FALSE)
  write.csv(go, "outputs/GO_results.csv")
}

# -------- Differential variability --------
if (isTRUE(steps$dv)) {
  logmsg("Differential variability (limma varFit)...")
  Mwork <- if (isTRUE(steps$filter)) readRDS("outputs/M_filtered.rds") else readRDS("outputs/M_raw.rds")
  g <- factor(targets[[cfg$model$group_column %||% "Sample_Group"]])
  design_dv <- model.matrix(~ g)
  fitvar <- varFit(Mwork, design=design_dv, coef=2)
  dv <- topVar(fitvar, number=Inf, adjust.method="BH")
  write.csv(dv, "outputs/DV_results.csv")
}

# -------- Cell-type composition note --------
if (isTRUE(steps$celltype_note)) {
  writeLines("Cell-type deconvolution requires a tissue-appropriate reference; standard Houseman (blood) is not suitable for solid tissues like kidney. Consider reference-free or organ-specific references.", "outputs/celltype_note.txt")
}

sink("outputs/session_info.txt"); print(sessionInfo()); sink()
logmsg("Done. Check outputs/.")
