# run_all.R — methyCleanr (universal; with sample selection + structured outputs)

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

# -------- OUTPUT LAYOUT --------
OUT <- list(
  qc   = "outputs/qc",
  data = "outputs/data",
  dmp  = "outputs/dmp",
  dmr  = "outputs/dmr",
  go   = "outputs/go",
  dv   = "outputs/dv",
  logs = "outputs/logs"
)
invisible(lapply(OUT, dir.create, recursive = TRUE, showWarnings = FALSE))

# Simple logger to outputs/logs/run_all.log
log_file_con <- file(file.path(OUT$logs, "run_all.log"), open = "wt")
sink(log_file_con, type = "output")
sink(log_file_con, type = "message")
.on_exit_log <- function() {
  try(sink(type="message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log_file_con), silent = TRUE)
}
on.exit(.on_exit_log(), add = TRUE)

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n")
save_rds <- function(obj, path) { saveRDS(obj, path); logmsg("Wrote %s", path) }
save_csv <- function(df, path, rn = FALSE) { utils::write.csv(df, path, row.names = rn); logmsg("Wrote %s", path) }
save_png <- function(expr, path, w=1200, h=700) { png(path, width=w, height=h); force(expr); dev.off(); logmsg("Wrote %s", path) }

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
  if (length(miss)) stop("Missing columns in samples.csv: ", paste(miss, collapse=", "))
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
    if (!idcol %in% names(pheno)) stop("ids_column '", idcol, "' not found in samples.csv.")
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
# default to samples.csv (your repo uses this naming)
pheno <- read.csv(input$phenotype_csv %||% "config/samples.csv",
                  stringsAsFactors = FALSE, check.names = FALSE)
pheno <- normalize_pheno(pheno)

# Apply selection BEFORE building targets
pheno <- apply_selection(pheno, select)
if (nrow(pheno) == 0) {
  stop("No rows selected from samples.csv. Check your select block in config/config.yaml.")
}

if (mode == "idat") {
  require_cols(pheno, c("Slide","Array"))
  if (!"Sample_Name" %in% names(pheno)) pheno$Sample_Name <- seq_len(nrow(pheno))
  if (!"Sample_Group" %in% names(pheno)) pheno$Sample_Group <- "Group"
  targets <- pheno[, intersect(c("Sample_Name","Sample_Group","Slide","Array"), names(pheno)), drop = FALSE]
  save_csv(targets, file.path(OUT$qc, "sample_sheet_used.csv"))
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
  save_rds(rgSet, file.path(OUT$data, "rgSet.rds"))

} else if (mode == "rgset") {
  rg_path <- input$rgset_rds %||% "config/rgSet.rds"
  if (!file.exists(rg_path)) stop("rgSet.rds not found at: ", rg_path)
  rgSet <- readRDS(rg_path)
  require_cols(pheno, c("Sample_Name"))
  if (!"Sample_Group" %in% names(pheno)) pheno$Sample_Group <- "Group"
  sm <- sampleNames(rgSet)
  rownames(pheno) <- pheno$Sample_Name
  if (!all(sm %in% rownames(pheno))) {
    stop("samples.csv must contain all sample names present in rgSet (Sample_Name column).")
  }
  pheno <- pheno[sm, , drop=FALSE]
  targets <- pheno[, intersect(c("Sample_Name","Sample_Group","Slide","Array"), names(pheno)), drop = FALSE]
  save_csv(targets, file.path(OUT$qc, "sample_sheet_used.csv"))
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
  save_rds(detP, file.path(OUT$data, "detectionP.rds"))

  save_png({
    barplot(colMeans(detP, na.rm = TRUE), las = 2, ylab = "Mean detection P-value")
    abline(h = qc$detP_threshold %||% 0.01, lty = 2)
  }, file.path(OUT$qc, "mean_detection_p_per_sample.png"))

  try(qcReport(rgSet, sampNames = targets$Sample_Name, sampGroups = targets$Sample_Group,
               pdf = file.path(OUT$qc, "qcReport.pdf")), silent = TRUE)
}

# Compute per-sample QC metrics
  thr <- qc$detP_threshold %||% 0.01
  sample_mean_detP <- colMeans(detP, na.rm = TRUE)
  sample_failed_frac <- colMeans(detP > thr, na.rm = TRUE)

# Decide samples to drop
  bad_mean   <- sample_mean_detP > (qc$sample_mean_detP_max %||% 0.01)
  bad_failed <- sample_failed_frac > (qc$sample_max_failed_frac %||% 0.05)

  drop_idx <- which(bad_mean | bad_failed)
  qc_tbl <- data.frame(
    Sample = colnames(detP),
    mean_detP = sample_mean_detP,
    failed_frac = sample_failed_frac,
    drop_mean = bad_mean,
    drop_failed = bad_failed,
    drop = (bad_mean | bad_failed)
  )

# Save QC summary + plot
  save_csv(qc_tbl, file.path(OUT$qc, "sample_qc_summary.csv"))
  save_png({
    op <- par(mfrow=c(1,2), mar=c(8,4,2,1))
    barplot(sample_mean_detP, las=2, ylab="Mean detP"); abline(h=thr, lty=2, col="red")
    barplot(sample_failed_frac, las=2, ylab=sprintf("Frac detP > %.2g", thr)); abline(h=(qc$sample_max_failed_frac %||% 0.05), lty=2, col="red")
    par(op)
  }, file.path(OUT$qc, "sample_qc_bars.png"))

# Apply sample filtering if needed
  if (length(drop_idx) > 0) {
    logmsg("Dropping %d low-quality samples: %s", length(drop_idx), paste(colnames(detP)[drop_idx], collapse=", "))
    keep_samples <- setdiff(seq_len(ncol(rgSet)), drop_idx)
    rgSet   <- rgSet[, keep_samples]
    detP    <- detP[, keep_samples, drop=FALSE]
    targets <- targets[keep_samples, , drop=FALSE]
  }

# -------- Normalization --------
if (isTRUE(steps$normalize)) {
  norm_method <- tolower(cfg$normalization$method %||% "noob_quantile")
  logmsg("Normalizing with method: %s ...", norm_method)

  if (norm_method == "funnorm") {
    mSet <- preprocessFunnorm(rgSet)
  } else {
    mSet.noob <- preprocessNoob(rgSet)
    mSet <- preprocessQuantile(mSet.noob)
  }

  save_rds(mSet, file.path(OUT$data, "mSet.rds"))
  beta <- getBeta(mSet); M <- getM(mSet)
  save_rds(beta, file.path(OUT$data, "beta_raw.rds"))
  save_rds(M,    file.path(OUT$data, "M_raw.rds"))
} else {
  mSet <- readRDS(file.path(OUT$data, "mSet.rds"))
  beta <- readRDS(file.path(OUT$data, "beta_raw.rds"))
  M    <- readRDS(file.path(OUT$data, "M_raw.rds"))
}

# -------- Explore --------
if (isTRUE(steps$explore)) {
  logmsg("Generating density & MDS plots...")
  save_png({ minfi::plotDensities(beta, main = "Beta densities (post-normalization)") },
           file.path(OUT$qc, "beta_density.png"))
  save_png({
    plotMDS(getM(mSet), top = 1000, gene.selection = "common",
            labels = targets$Sample_Name, col = as.numeric(factor(targets$Sample_Group)))
  }, file.path(OUT$qc, "mds.png"))
}

# -------- Filtering --------
if (isTRUE(steps$filter)) {
  logmsg("Filtering probes...")
  detP <- if (exists("detP")) detP else readRDS(file.path(OUT$data, "detectionP.rds"))
  keep_det   <- rowMeans(detP < (qc$detP_threshold %||% 0.01), na.rm = TRUE) >= (qc$detP_fraction_pass %||% 0.90)
  keep_nonNA <- rowMeans(is.na(beta)) <= (qc$max_na_fraction %||% 0.05)

  ann <- getAnnotation(ann_obj)
  sex_mask <- ann$chr %in% c("chrX","chrY")
  keep_sex <- if (isTRUE(qc$remove_sex)) !sex_mask else rep(TRUE, length(sex_mask))

  idx <- match(rownames(beta), rownames(ann))
  keep <- keep_det & keep_nonNA & keep_sex[idx]

  if ((extra$cross_reactive_list %||% "") != "") {
    xr <- unique(trimws(readLines(extra$cross_reactive_list)))
    keep <- keep & !(rownames(beta) %in% xr)
  }
  if ((extra$snp_probes_list %||% "") != "") {
    sp <- unique(trimws(readLines(extra$snp_probes_list)))
    keep <- keep & !(rownames(beta) %in% sp)
  }

  beta_f <- beta[keep,,drop=FALSE]; M_f <- M[keep,,drop=FALSE]
  save_rds(beta_f, file.path(OUT$data, "beta_filtered.rds"))
  save_rds(M_f,    file.path(OUT$data, "M_filtered.rds"))

  utils::write.table(rownames(beta_f), file.path(OUT$data, "filtered_probe_ids.txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  save_csv(data.frame(probes_kept = sum(keep)), file.path(OUT$qc, "filter_summary.csv"))
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
  Mwork <- if (isTRUE(steps$filter)) readRDS(file.path(OUT$data, "M_filtered.rds")) else readRDS(file.path(OUT$data, "M_raw.rds"))
  dd <- get_design(targets, cfg$model$group_column %||% "Sample_Group")
  cont <- parse_contrast(cfg$model$contrast %||% "Case- Control", dd$design)

  fit  <- lmFit(Mwork, dd$design)
  fit2 <- contrasts.fit(fit, cont); fit2 <- eBayes(fit2)
  dmp  <- topTable(fit2, number = Inf, adjust.method = "BH")
  utils::write.csv(dmp, file.path(OUT$dmp, "DMP_results.csv"), row.names = TRUE)
  logmsg("Wrote %s", file.path(OUT$dmp, "DMP_results.csv"))

  save_png({
    with(dmp, {
      plot(logFC, -log10(P.Value), pch=20, xlab="logFC (M)", ylab="-log10 P")
      abline(v=c(-1,1), lty=2); abline(h=1.3, lty=2)
    })
  }, file.path(OUT$dmp, "volcano_dmp.png"), w=1100, h=800)
}

# -------- DMR --------
if (isTRUE(steps$dmr)) {
  logmsg("Regional DMRs (DMRcate)...")
  Mwork <- if (isTRUE(steps$filter)) readRDS(file.path(OUT$data, "M_filtered.rds")) else readRDS(file.path(OUT$data, "M_raw.rds"))
  dd <- get_design(targets, cfg$model$group_column %||% "Sample_Group")
  cont <- parse_contrast(cfg$model$contrast %||% "Case- Control", dd$design)
  annDf <- getAnnotation(ann_obj)

  myAnnot <- cpg.annotate(object="array", datatype="array", what="M", analysis.type="differential",
                          x=Mwork, design=dd$design, contrasts=TRUE, cont.matrix=cont, coef=1,
                          fdr=0.05, annotation=annDf)
  dmrcoutput <- dmrcate(myAnnot, lambda=1000, C=2)
  ranges <- extractRanges(dmrcoutput, genome="hg19")
  save_rds(ranges, file.path(OUT$dmr, "DMR_ranges.rds"))

  dmr_df <- tryCatch({
    data.frame(
      chr     = as.character(GenomicRanges::seqnames(ranges)),
      start   = GenomicRanges::start(ranges),
      end     = GenomicRanges::end(ranges),
      width   = GenomicRanges::width(ranges),
      no_cpgs = GenomicRanges::mcols(ranges)[,"no.cpgs"],
      mean_fc = suppressWarnings(GenomicRanges::mcols(ranges)[,"meanbetafc"]),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
  if (!is.null(dmr_df)) save_csv(dmr_df, file.path(OUT$dmr, "DMR_summary.csv"))
}

# -------- GO --------
if (isTRUE(steps$go)) {
  logmsg("GO testing (gometh)...")
  dmp <- utils::read.csv(file.path(OUT$dmp, "DMP_results.csv"), stringsAsFactors = FALSE)
  Mwork <- if (isTRUE(steps$filter)) readRDS(file.path(OUT$data, "M_filtered.rds")) else readRDS(file.path(OUT$data, "M_raw.rds"))
  sig <- dmp$X[dmp$adj.P.Val < 0.05]
  go <- gometh(sig.cpg = sig, all.cpg = rownames(Mwork), collection="GO",
               array.type = ifelse(ann_choice=='EPIC','EPIC','450K'), plot.bias=FALSE)
  utils::write.csv(go, file.path(OUT$go, "GO_results.csv"), row.names = TRUE)
  logmsg("Wrote %s", file.path(OUT$go, "GO_results.csv"))
}

# -------- Differential variability --------
if (isTRUE(steps$dv)) {
  logmsg("Differential variability (limma varFit)...")
  Mwork <- if (isTRUE(steps$filter)) readRDS(file.path(OUT$data, "M_filtered.rds")) else readRDS(file.path(OUT$data, "M_raw.rds"))
  g <- factor(targets[[cfg$model$group_column %||% "Sample_Group"]])
  design_dv <- model.matrix(~ g)
  fitvar <- varFit(Mwork, design=design_dv, coef=2)
  dv <- topVar(fitvar, number=Inf, adjust.method="BH")
  utils::write.csv(dv, file.path(OUT$dv, "DV_results.csv"), row.names = TRUE)
  logmsg("Wrote %s", file.path(OUT$dv, "DV_results.csv"))
}

# -------- Cell-type composition note --------
if (isTRUE(steps$celltype_note)) {
  writeLines("Cell-type deconvolution requires a tissue-appropriate reference; standard Houseman (blood) is not suitable for solid tissues like kidney. Consider reference-free or organ-specific references.",
             file.path(OUT$logs, "celltype_note.txt"))
}

# -------- Session info --------
utils::writeLines(capture.output(sessionInfo()), con = file.path(OUT$logs, "session_info.txt"))
logmsg("Done. Check outputs/.")
