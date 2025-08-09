pkgs <- c(
  "limma",
  "minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kmanifest",
  "RColorBrewer",
  "missMethyl",
  "matrixStats",
  "minfiData",
  "Gviz",
  "DMRcate",
  "stringr"
)

BiocManager::install(pkgs, ask = FALSE, update = FALSE)
