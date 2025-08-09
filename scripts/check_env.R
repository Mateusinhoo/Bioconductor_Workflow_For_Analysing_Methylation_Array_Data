pkgs <- c("minfi","missMethyl","DMRcate","Gviz",
          "IlluminaHumanMethylation450kmanifest",
          "IlluminaHumanMethylation450kanno.ilmn12.hg19",
          "BiocParallel","Rsamtools","limma")

for (p in pkgs) {
  cat(sprintf("%-45s : %s\n", p, tryCatch(as.character(packageVersion(p)),
                                          error=function(e) "NOT INSTALLED")))
}

if (requireNamespace("sessioninfo", quietly = TRUE)) {
  sessioninfo::session_info(pkgs = pkgs)
} else {
  cat("\nTip: install.packages('sessioninfo') for a full dump.\n")
}
