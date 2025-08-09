# scripts/sanitize_idats.R
# Usage: source("scripts/sanitize_idats.R"); sanitize_idats("data/idats_raw","data/idats_clean")

sanitize_idats <- function(in_dir, out_dir = file.path(in_dir, "..", "idats_clean")) {
  if (!dir.exists(in_dir)) stop("Input directory not found: ", in_dir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  files <- list.files(in_dir, pattern = "(?i)idat(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop("No IDAT files found under: ", in_dir)

  n_copied <- 0L
  for (f in files) {
    bn <- basename(f)

    # Normalize extension to lowercase .idat[.gz]
    bn <- sub("\\.IDAT(\\.GZ)?$", ".idat\\1", bn, ignore.case = TRUE)

    # Normalize channel hints _R/_G to _Red/_Grn
    if (grepl("_R\\.", bn, ignore.case = TRUE)) bn <- sub("_R\\.", "_Red.", bn)
    if (grepl("_G\\.", bn, ignore.case = TRUE)) bn <- sub("_G\\.", "_Grn.", bn)

    dest <- file.path(out_dir, bn)
    if (!file.exists(dest)) {
      ok <- file.copy(f, dest, overwrite = FALSE)
      if (ok) n_copied <- n_copied + 1L
    }
  }
  message("Sanitized copies written to: ", out_dir, " (", n_copied, " files)")
  return(invisible(out_dir))
}
