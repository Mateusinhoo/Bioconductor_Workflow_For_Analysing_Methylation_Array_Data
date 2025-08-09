# scripts/check_env.R
# Quick diagnostics to confirm your R + package environment

cat("==== R sessionInfo ====\n")
print(sessionInfo())

cat("\n==== Key packages ====\n")
pkgs <- c(
  "yaml","minfi","limma","DMRcate","missMethyl",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "ggplot2","matrixStats","R.utils","optparse"
)
for (p in pkgs) {
  cat(sprintf("%-45s : %s\n", p,
              if (suppressWarnings(requireNamespace(p, quietly = TRUE)))
                as.character(utils::packageVersion(p)) else "NOT INSTALLED"))
}

# Optional: quick GPU/BLAS hints
cat("\n==== BLAS/LAPACK vendor (if shown) ====\n")
cap <- capabilities()
print(cap)
