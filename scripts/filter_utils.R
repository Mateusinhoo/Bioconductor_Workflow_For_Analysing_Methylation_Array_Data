# scripts/filter_utils.R
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

get_annotation_for <- function(ann_choice = c("450K","EPIC")) {
  ann_choice <- match.arg(ann_choice)
  if (ann_choice == "EPIC") {
    getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  } else {
    getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
}

mask_detection <- function(detP, threshold = 0.01, fraction_pass = 0.90) {
  rowMeans(detP < threshold, na.rm = TRUE) >= fraction_pass
}

mask_nonNA <- function(beta, max_na_fraction = 0.05) {
  rowMeans(is.na(beta)) <= max_na_fraction
}

mask_sex_chrom <- function(ann_df, drop_sex = FALSE) {
  if (!drop_sex) return(rep(TRUE, nrow(ann_df)))
  !(ann_df$chr %in% c("chrX", "chrY"))
}

apply_probe_lists <- function(keep, beta_rownames, cross_reactive_file = "", snp_file = "") {
  if (nzchar(cross_reactive_file) && file.exists(cross_reactive_file)) {
    xr <- unique(trimws(readLines(cross_reactive_file)))
    keep <- keep & !(beta_rownames %in% xr)
  }
  if (nzchar(snp_file) && file.exists(snp_file)) {
    sp <- unique(trimws(readLines(snp_file)))
    keep <- keep & !(beta_rownames %in% sp)
  }
  keep
}
