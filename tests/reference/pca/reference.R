# What R's prcomp gives for the matrices that make_reference.py wrote.
# To standardize, prcomp divides by a standard deviation that has n - 1 in
# it and pyNei by one that has n, so the projections of a standardized
# prcomp are multiplied by sqrt(n / (n - 1)) to be pyNei's. Without
# standardizing the two give the same projections.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
cat(R.version.string, "\n")

fix_signs <- function(p) {
  for (k in seq_len(ncol(p$x))) {
    s <- sign(p$x[which.max(abs(p$x[, k])), k])
    p$x[, k] <- p$x[, k] * s
    p$rotation[, k] <- p$rotation[, k] * s
  }
  p
}

write_pca <- function(name, x, scale, num_comps) {
  n <- nrow(x)
  p <- fix_signs(prcomp(x, center = TRUE, scale. = scale))
  k <- min(num_comps, ncol(p$x))
  write.table(format(p$x[, 1:k, drop = FALSE] * (if (scale) sqrt(n / (n - 1)) else 1), digits = 12), file.path(out, paste0(name, ".r.projections.tsv")), sep = "\t", quote = FALSE)
  write.table(format((p$sdev^2 / sum(p$sdev^2) * 100)[1:k], digits = 12), file.path(out, paste0(name, ".r.percent.tsv")), sep = "\t", quote = FALSE, col.names = FALSE)
  write.table(format(t(p$rotation[, 1:k, drop = FALSE]), digits = 12), file.path(out, paste0(name, ".r.princomps.tsv")), sep = "\t", quote = FALSE)
}

dosages <- function(name) {
  m <- as.matrix(read.table(file.path(out, paste0(name, ".mat012.tsv")), sep = "\t"))
  m[m == -1] <- NA
  x <- t(m)
  colnames(x) <- seq_len(ncol(x)) - 1
  for (j in seq_len(ncol(x))) x[is.na(x[, j]), j] <- mean(x[, j], na.rm = TRUE)
  x[, apply(x, 2, sd) > 0, drop = FALSE]
}

write_pca("sim_missing", dosages("sim_missing"), TRUE, 10)
write_pca("worked", dosages("worked"), TRUE, 4)
write_pca("worked3", dosages("worked3"), TRUE, 4)
iris4 <- as.matrix(read.table(file.path(out, "iris.tsv"), sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
write_pca("iris", iris4, TRUE, 4)
write_pca("iris_not_standardized", iris4, FALSE, 4)
