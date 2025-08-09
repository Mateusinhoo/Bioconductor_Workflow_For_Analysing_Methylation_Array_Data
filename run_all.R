# run_all.R — Kidney 450k workflow (KIRC)
# Toggle to do cleaning only (up through filtering), or full EWAS
CLEAN_ONLY <- TRUE

# ---------------- packages & project setup ----------------
message("Loading packages...")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c(
  "minfi","limma","IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19","DMRcate","missMethyl",
  "matrixStats","ggplot2","R.utils"
)
for (p in pkgs) {
  if (!suppressWarnings(requireNamespace(p, quietly = TRUE))) {
    if (p %in% rownames(installed.packages())) next
    if (p %in% c("minfi","DMRcate","missMethyl","IlluminaHumanMethylation450kmanifest","IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
      BiocManager::install(p, ask = FALSE, update = FALSE)
    } else {
      install.packages(p, dependencies = TRUE)
    }
  }
}
suppressPackageStartupMessages({
  library(minfi); library(limma); library(DMRcate); library(missMethyl)
  library(matrixStats); library(ggplot2); library(R.utils)
})

dir.create("outputs/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/normalized", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/results", recursive = TRUE, showWarnings = FALSE)

# ---------------- user inputs ----------------
idat_dir <- "data/idats_raw"                        # put your IDATs here
pheno_csv <- "config/phenotype_kirc.csv"            # provided (443 KIRC rows)
sample_sheet_out <- "config/sample_sheet.csv"       # derived for minfi
project_name <- "KIRC"

# ---------------- SampleSheet derivation ----------------
stopifnot(file.exists(pheno_csv))
pheno <- read.csv(pheno_csv, stringsAsFactors = FALSE, check.names = FALSE)

# Keep only KIRC rows (in case file expands later)
if ("Study" %in% names(pheno)) {
  pheno <- subset(pheno, grepl("Kidney renal clear cell carcinoma", Study, fixed = TRUE))
}

# Require the Illumina IDs & positions
id_cols <- c("Sentrix ID","Sentrix position")
missing_cols <- setdiff(id_cols, names(pheno))
if (length(missing_cols)) stop("Missing columns in phenotype CSV: ", paste(missing_cols, collapse=", "))

# Build a minimal SampleSheet understood by minfi
# minfi likes columns Slide (Sentrix_ID) and Array (SentrixPosition)
Sample_Name <- if ("TCGA ID" %in% names(pheno)) pheno[["TCGA ID"]] else seq_len(nrow(pheno))
Sample_Group <- if ("Type" %in% names(pheno)) pheno[["Type"]] else "Unknown"
Slide <- pheno[["Sentrix ID"]]
Array <- pheno[["Sentrix position"]]

sample_sheet <- data.frame(
  Sample_Name = Sample_Name,
  Sample_Group = Sample_Group,
  Slide = Slide,
  Array = Array,
  stringsAsFactors = FALSE
)

write.csv(sample_sheet, sample_sheet_out, row.names = FALSE)
message("Wrote derived SampleSheet to: ", sample_sheet_out)

# ---------------- Load IDATs ----------------
message("Reading IDATs... (pointing at idat_dir = ", idat_dir, ")")
if (!dir.exists(idat_dir)) stop("IDAT directory not found: ", idat_dir)

# read.metharray.sheet expects the SampleSheet in the base path; we pass explicit targets instead.
targets <- read.csv(sample_sheet_out, stringsAsFactors = FALSE)
# Attach path hints for minfi: it will search recursively under idat_dir by Slide/Array naming
# If your files are non-standard named, set force=TRUE to let minfi scan recursively.
rgSet <- read.metharray.exp(base = idat_dir, targets = targets, extended = TRUE, force = TRUE)
saveRDS(rgSet, file = "outputs/normalized/rgSet.rds")

# ---------------- Quality control ----------------
message("Computing detection P-values...")
detP <- detectionP(rgSet)
saveRDS(detP, "outputs/qc/detectionP_matrix.rds")

# Mean detection P per sample
detP_means <- colMeans(detP, na.rm = TRUE)
png("outputs/qc/mean_detection_p_per_sample.png", width = 1200, height = 700)
barplot(detP_means, las = 2, ylab = "Mean detection P-value")
abline(h = 0.01, lty = 2)
dev.off()

# qcReport (tries to use targets$Sample_Group if available)
targets$ID <- targets$Sample_Name
try({
  qcReport(rgSet, sampNames = targets$ID, sampGroups = targets$Sample_Group,
           pdf = "outputs/qc/qcReport.pdf")
}, silent = TRUE)

# ---------------- Normalization ----------------
# Robust baseline: background correction + dye-bias (Noob), then between-array normalization (Quantile)
message("Normalizing with preprocessNoob -> preprocessQuantile ...")
mSet.noob <- preprocessNoob(rgSet)
mSet <- preprocessQuantile(mSet.noob)
saveRDS(mSet, file = "outputs/normalized/mSet_quantile_noob.rds")

beta <- getBeta(mSet)
M <- getM(mSet)
saveRDS(beta, "outputs/normalized/beta.rds")
saveRDS(M, "outputs/normalized/M.rds")

# QC density plots (beta)
png("outputs/qc/beta_density.png", width = 1200, height = 700)
plotDensities(beta, main = "Beta densities (post-normalization)")
dev.off()

# MDS for exploration
png("outputs/qc/mds_by_sample_group.png", width = 1200, height = 700)
plotMDS(M, top = 1000, gene.selection = "common", labels = targets$Sample_Name, col = as.numeric(factor(targets$Sample_Group)))
dev.off()

# ---------------- Filtering ----------------
# 1) Drop poor-quality probes (detection P > 0.01 in >10% samples)
keep_det <- rowMeans(detP < 0.01, na.rm = TRUE) > 0.90
# 2) Remove probes with NA across many samples
keep_nonNA <- rowMeans(is.na(beta)) < 0.05
# 3) Optionally remove sex-chr probes (toggle)
REMOVE_SEX <- FALSE
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
is_sex <- ann450k@data$chr %in% c("chrX","chrY")
keep_sex <- if (REMOVE_SEX) !is_sex else rep(TRUE, length(is_sex))

keep <- keep_det & keep_nonNA & keep_sex[match(rownames(beta), rownames(ann450k@data))]
beta_f <- beta[keep, , drop = FALSE]
M_f <- M[keep, , drop = FALSE]
saveRDS(keep, "outputs/normalized/filter_keep_mask.rds")
saveRDS(beta_f, "outputs/normalized/beta_filtered.rds")
saveRDS(M_f, "outputs/normalized/M_filtered.rds")

png("outputs/qc/num_probes_kept.png", width = 900, height = 600)
plot(c(sum(keep)), type="h", xaxt="n", xlab="", ylab="Probes kept"); axis(1, at=1, labels=FALSE)
dev.off()

if (CLEAN_ONLY) {
  message("CLEAN_ONLY=TRUE complete. Cleaned objects saved in outputs/normalized/.")
  quit(save = "no")
}

# ---------------- Design & contrasts (Tumor vs Normal) ----------------
group <- factor(targets$Sample_Group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
saveRDS(design, "outputs/models/design_matrix.rds")

# ---------------- Probe-wise DMP (limma) ----------------
fit <- lmFit(M_f, design)
cont <- makeContrasts(Tumor_vs_Normal = `Primary solid Tumor` - `Solid Tissue Normal`, levels = design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)
dmp <- topTable(fit2, number = Inf, adjust.method = "BH")
write.csv(dmp, "outputs/results/DMP_Tumor_vs_Normal.csv", row.names = TRUE)

# Volcano
png("outputs/results/volcano_dmp.png", width = 1100, height = 800)
with(dmp, {
  plot(logFC, -log10(P.Value), pch=20, xlab="logFC (M-scale)", ylab="-log10 P")
  abline(v = c(-1,1), lty = 2); abline(h = 1.3, lty = 2) # ~P=0.05
})
dev.off()

# ---------------- Regional DMRs (DMRcate) ----------------
myAnnot <- cpg.annotate(
  object = "array", datatype = "array", what = "M", analysis.type = "differential",
  x = M_f, design = design, contrasts = TRUE, cont.matrix = cont, coef = 1,
  fdr = 0.05, annotation = ann450k
)
dmrcoutput <- dmrcate(myAnnot, lambda = 1000, C = 2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
saveRDS(results.ranges, "outputs/results/DMR_ranges.rds")

# ---------------- GO testing (missMethyl::gometh) ----------------
sig <- rownames(dmp)[which(dmp$adj.P.Val < 0.05)]
go <- gometh(sig.cpg = sig, all.cpg = rownames(M_f), collection = "GO", array.type = "450K", plot.bias = FALSE)
write.csv(go, "outputs/results/GO_gometh.csv", row.names = TRUE)

# ---------------- Differential variability (limma varFit) ----------------
# Note: DV in tumors may highlight heterogeneity regions
design_dv <- model.matrix(~ group)
fitvar <- varFit(M_f, design = design_dv, coef = 2)
dv <- topVar(fitvar, number = Inf, adjust.method = "BH")
write.csv(dv, "outputs/results/DiffVariability_Tumor_vs_Normal.csv", row.names = TRUE)

# ---------------- Cell type composition (discussion) ----------------
# Placeholder: Houseman estimates are not appropriate for solid tumors.
writeLines("Reference-based cell-type deconvolution (Houseman) is blood-specific and not applicable to kidney tumor tissue. Consider reference-free methods or single-cell references if available.", 
           con = "outputs/results/cell_type_note.txt")

message("Full analysis complete. Results written to outputs/.")
