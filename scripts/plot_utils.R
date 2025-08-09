# scripts/plot_utils.R
suppressPackageStartupMessages({ library(ggplot2) })

plot_mean_detP <- function(detP, out_png, threshold = 0.01) {
  stopifnot(is.matrix(detP))
  m <- colMeans(detP, na.rm = TRUE)
  png(out_png, width = 1200, height = 700)
  barplot(m, las = 2, ylab = "Mean detection P-value")
  abline(h = threshold, lty = 2)
  dev.off()
  invisible(out_png)
}

plot_beta_density <- function(beta, out_png, title = "Beta densities") {
  png(out_png, width = 1200, height = 700)
  minfi::plotDensities(beta, main = title)
  dev.off()
  invisible(out_png)
}

plot_mds <- function(M, sample_names, groups, out_png) {
  png(out_png, width = 1200, height = 700)
  plotMDS(M, top = 1000, gene.selection = "common",
          labels = sample_names, col = as.numeric(factor(groups)))
  dev.off()
  invisible(out_png)
}

plot_volcano_limma <- function(tbl, out_png, vlines = c(-1, 1), h_p = 0.05) {
  png(out_png, width = 1100, height = 800)
  with(tbl, {
    plot(logFC, -log10(P.Value), pch = 20, xlab = "logFC (M)", ylab = "-log10 P")
    abline(v = vlines, lty = 2)
    abline(h = -log10(h_p), lty = 2)
  })
  dev.off()
  invisible(out_png)
}
