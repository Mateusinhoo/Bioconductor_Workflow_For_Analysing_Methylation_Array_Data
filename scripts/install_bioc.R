light_pkgs <- c(
  "RColorBrewer",  # plotting
  "stringr"        # already in env

if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Prefer macOS binaries from CRAN when possible
if (Sys.info()[["sysname"]] == "Darwin") {
  options(pkgType = "mac.binary")
}

need <- setdiff(light_pkgs, rownames(installed.packages()))
if (length(need)) {
  BiocManager::install(need, ask = FALSE, update = FALSE,
                       Ncpus = max(1, parallel::detectCores() - 1))
} else {
  message("All requested CRAN/Bioc light packages already installed.")
}
