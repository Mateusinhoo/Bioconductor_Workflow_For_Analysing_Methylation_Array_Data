# run_all.R — methyCleanr (universal)
suppressPackageStartupMessages({
  library(yaml); library(minfi); library(limma); library(DMRcate); library(missMethyl)
  library(ggplot2); library(matrixStats); library(R.utils)
})

cfg <- yaml::read_yaml("config/config.yaml")
prj  <- cfg$project_name
idat_dir <- cfg$idat_dir
pheno_csv <- cfg$phenotype_csv
steps <- cfg$steps
qc <- cfg$qc
array_type <- tolower(cfg$array_type %||% "auto")

`%||%` <- function(a,b) if (!is.null(a)) a else b

dir.create("outputs", showWarnings = FALSE)

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n")

# ---- Read phenotype / build targets ----
stopifnot(file.exists(pheno_csv))
pheno <- read.csv(pheno_csv, stringsAsFactors = FALSE, check.names = FALSE)

# normalize column names
n <- names(pheno)
n <- sub("^Sentrix ID$","Slide", n, ignore.case = TRUE)
n <- sub("^Sentrix position$","Array", n, ignore.case = TRUE)
names(pheno) <- n

if (!all(c("Slide","Array") %in% names(pheno))) {
  stop("Phenotype CSV must include Slide and Array columns (Sentrix ID and Sentrix position).")
}
if (!"Sample_Name" %in% names(pheno)) pheno$Sample_Name <- seq_len(nrow(pheno))
if (!"Sample_Group" %in% names(pheno)) pheno$Sample_Group <- "Group"

targets <- pheno[, c("Sample_Name","Sample_Group","Slide","Array")]
write.csv(targets, "outputs/sample_sheet_used.csv", row.names = FALSE)

# ---- Load ----
if (isTRUE(steps$load)) {
  logmsg("Loading IDATs from %s ...", idat_dir)
  if (!dir.exists(idat_dir)) stop("IDAT directory does not exist: ", idat_dir)
  rgSet <- read.metharray.exp(base = idat_dir, targets = targets, extended = TRUE, force = TRUE)
  saveRDS(rgSet, file = "outputs/rgSet.rds")
  logmsg("rgSet saved.")
} else {
  rgSet <- readRDS("outputs/rgSet.rds")
}

# ---- Array type detection (if auto) ----
ann_choice <- NULL
if (array_type == "auto") {
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

# ---- QC ----
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

# ---- Normalization ----
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

# ---- Explore ----
if (isTRUE(steps$explore)) {
  logmsg("Making density/MDS plots...")
  png("outputs/beta_density.png", width = 1200, height = 700)
  plotDensities(beta, main = "Beta densities (post-normalization)")
  dev.off()
  png("outputs/mds.png", width = 1200, height = 700)
  plotMDS(getM(mSet), top = 1000, gene.selection = "common",
          labels = targets$Sample_Name, col = as.numeric(factor(targets$Sample_Group)))
  dev.off()
}

# ---- Filter ----
if (isTRUE(steps$filter)) {
  logmsg("Filtering probes...")
  detP <- if (exists("detP")) detP else readRDS("outputs/detectionP.rds")
  keep_det <- rowMeans(detP < (qc$detP_threshold %||% 0.01), na.rm = TRUE) >= (qc$detP_fraction_pass %||% 0.90)
  keep_nonNA <- rowMeans(is.na(beta)) <= (qc$max_na_fraction %||% 0.05)
  ann <- getAnnotation(ann_obj)
  sex_mask <- ann$chr %in% c("chrX","chrY")
  keep_sex <- if (isTRUE(qc$remove_sex)) !sex_mask else rep(TRUE, length(sex_mask))
  # align masks
  idx <- match(rownames(beta), rownames(ann))
  keep <- keep_det & keep_nonNA & keep_sex[idx]
  beta <- beta[keep,,drop=FALSE]; M <- M[keep,,drop=FALSE]
  saveRDS(beta, "outputs/beta_filtered.rds")
  saveRDS(M, "outputs/M_filtered.rds")
  write.csv(data.frame(probes_kept=sum(keep)), "outputs/filter_summary.csv", row.names = FALSE)
}

# ---- Modeling helpers ----
get_design <- function(targets, group_col) {
  g <- factor(targets[[group_col]])
  design <- model.matrix(~ 0 + g); colnames(design) <- levels(g)
  list(g=g, design=design)
}
parse_contrast <- function(contrast, design) {
  makeContrasts(contrasts = contrast, levels = design)
}

# ---- DMP ----
if (isTRUE(steps$dmp)) {
  logmsg("Probe-wise DMP with limma...")
  Mfit <- if (isTRUE(steps$filter)) readRDS("outputs/M_filtered.rds") else readRDS("outputs/M_raw.rds")
  dd <- get_design(targets, cfg$model$group_column %||% "Sample_Group")
  cont <- parse_contrast(cfg$model$contrast %||% "Case- Control", dd$design)
  fit <- lmFit(Mfit, dd$design)
  fit2 <- contrasts.fit(fit, cont); fit2 <- eBayes(fit2)
  dmp <- topTable(fit2, number = Inf, adjust.method = "BH")
  write.csv(dmp, "outputs/DMP_results.csv")
  png("outputs/volcano_dmp.png", width = 1100, height = 800)
  with(dmp, { plot(logFC, -log10(P.Value), pch=20, xlab="logFC (M)", ylab="-log10 P"); abline(v=c(-1,1), lty=2); abline(h=1.3, lty=2) })
  dev.off()
}

# ---- DMR ----
if (isTRUE(steps$dmr)) {
  logmsg("Regional DMRs with DMRcate...")
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

# ---- GO ----
if (isTRUE(steps$go)) {
  logmsg("GO testing with gometh...")
  dmp <- read.csv("outputs/DMP_results.csv", stringsAsFactors = FALSE)
  Mwork <- if (isTRUE(steps$filter)) readRDS("outputs/M_filtered.rds") else readRDS("outputs/M_raw.rds")
  sig <- rownames(Mwork)[rownames(Mwork) %in% dmp$X & dmp$adj.P.Val < 0.05]
  go <- gometh(sig.cpg = sig, all.cpg = rownames(Mwork), collection="GO",
               array.type = ifelse(ann_choice=='EPIC','EPIC','450K'), plot.bias=FALSE)
  write.csv(go, "outputs/GO_results.csv")
}

# ---- Differential variability ----
if (isTRUE(steps$dv)) {
  logmsg("Differential variability with limma varFit...")
  Mwork <- if (isTRUE(steps$filter)) readRDS("outputs/M_filtered.rds") else readRDS("outputs/M_raw.rds")
  g <- factor(targets[[cfg$model$group_column %||% "Sample_Group"]])
  design_dv <- model.matrix(~ g)
  fitvar <- varFit(Mwork, design=design_dv, coef=2)
  dv <- topVar(fitvar, number=Inf, adjust.method="BH")
  write.csv(dv, "outputs/DV_results.csv")
}

# ---- Cell type composition note ----
if (isTRUE(steps$celltype)) {
  writeLines("Cell-type deconvolution requires a tissue-appropriate reference; standard Houseman (blood) is not suitable for solid tissues like kidney. Consider reference-free or organ-specific references.", "outputs/celltype_note.txt")
}

logmsg("Done. Check outputs/.")
