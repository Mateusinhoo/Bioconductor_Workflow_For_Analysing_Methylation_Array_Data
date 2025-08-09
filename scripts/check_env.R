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

cat("R:", R.version.string, "\n")
ok <- sapply(pkgs, function(p) suppressPackageStartupMessages(require(p, character.only = TRUE, quietly = TRUE)))
print(data.frame(package = pkgs, loaded = ok), row.names = FALSE)

if (!all(ok)) quit(status = 1)
cat("Environment looks good.\n")
