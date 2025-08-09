# scripts/validate_samples.R
# Usage (inside R): source("scripts/validate_samples.R"); validate_samples("config/samples.csv")

validate_samples <- function(csv_path, required = c("Slide","Array"),
                             recommended = c("Sample_Name","Sample_Group")) {
  if (!file.exists(csv_path)) stop("File not found: ", csv_path)
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Map Sentrix ID/position -> Slide/Array if needed
  nm <- names(df)
  nm <- sub("^Sentrix ID$","Slide", nm, ignore.case = TRUE)
  nm <- sub("^Sentrix position$","Array", nm, ignore.case = TRUE)
  names(df) <- nm

  miss <- setdiff(required, names(df))
  if (length(miss)) stop("Missing required column(s) in samples.csv: ", paste(miss, collapse = ", "))

  # Basic checks
  if (anyNA(df$Slide) || anyNA(df$Array)) stop("Slide/Array contain NA values.")
  if (!all(grepl("^[0-9]+$", df$Slide))) warning("Some Slide values are not numeric-like.")
  if (!all(grepl("^[RC][0-9]{2}C[0-9]{2}$", df$Array))) warning("Some Array values don't match typical RxxCxx pattern.")

  # Hints
  rec_miss <- setdiff(recommended, names(df))
  if (length(rec_miss)) message("Note: missing recommended column(s): ", paste(rec_miss, collapse = ", "))

  message("samples.csv looks OK: ", nrow(df), " rows, ", ncol(df), " columns.")
  invisible(df)
}
