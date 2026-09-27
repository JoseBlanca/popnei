# What ape's pcoa gives for three distance matrices, the literals of the
# principal coordinates of docs/specs/pca.md. Run from the root of the
# repository with R 4.6.1 and ape 5.8.1:
#
#     Rscript tests/reference/pca/pcoa_reference.R tests/reference/pca
#
# The three matrices:
#
# - small, the ten distances of five individuals of test_pcoa of pyNei's
#   test/test_pca.py, which are not Euclidean: one eigenvalue is negative.
# - panel, the Kosman distances of the 200 individuals of
#   tests/reference/dists/panel.vcf.gz, as R's gd.kosman gave them in
#   tests/reference/dists/panel.gdkosman.tsv.
# - four_alleles, the Kosman distances of the 40 individuals of
#   tests/reference/dists/four_alleles.vcf.gz, from four_alleles.gdkosman.tsv,
#   which have no negative eigenvalue.
#
# pcoa with no correction gives the axes of the positive eigenvalues, and
# its Relative_eig is each eigenvalue over the sum of all of them, negative
# ones included. An eigenvalue counts as positive or negative when it is
# further from 0 than the largest times n times the machine epsilon, the
# threshold of the spec. Each axis gets the sign of the spec: its coordinate
# of the largest absolute value is positive. No axis of these matrices has
# two coordinates within 64 units in the last place of each other at its
# largest absolute value, which the script checks.
#
# It writes, for each matrix <name>:
#
# - <name>.pcoa.r.projections.tsv, one row per individual, one column per
#   axis;
# - <name>.pcoa.r.percent.tsv, the percentage of each axis;
# - <name>.pcoa.r.negative.tsv, the negative eigenvalues percent, 100 times
#   the sum of the absolute values of the negative eigenvalues over the sum
#   of all of them.
suppressMessages(library(ape))
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
if (as.character(getRversion()) != "4.6.1") stop("this script needs R 4.6.1")
if (as.character(packageVersion("ape")) != "5.8.1") stop("this script needs ape 5.8.1")
cat(R.version.string, ", ape ", as.character(packageVersion("ape")), "\n", sep = "")

square_of <- function(dist_vector) {
  n <- (1 + sqrt(1 + 8 * length(dist_vector))) / 2
  d <- matrix(0, n, n)
  # lower.tri walks down each column in turn, so it takes the pairs in the
  # order (0, 1), (0, 2), ..., (1, 2), ..., the order of dist_vector.
  d[lower.tri(d)] <- dist_vector
  d + t(d)
}

write_pcoa <- function(name, dist_vector, names) {
  d <- square_of(dist_vector)
  n <- nrow(d)
  p <- pcoa(as.dist(d), correction = "none")
  eig <- p$values$Eigenvalues
  tolerance <- max(eig) * n * .Machine$double.eps
  positive <- eig > tolerance
  negative <- eig < -tolerance
  coords <- p$vectors[, seq_len(sum(positive)), drop = FALSE]
  for (k in seq_len(ncol(coords))) {
    sorted <- sort(abs(coords[, k]), decreasing = TRUE)
    if (length(sorted) > 1 && sorted[1] - sorted[2] <= 64 * .Machine$double.eps * sorted[1]) {
      stop(sprintf("axis %d of %s has a tie at its largest absolute value", k, name))
    }
    coords[, k] <- coords[, k] * sign(coords[which.max(abs(coords[, k])), k])
  }
  rownames(coords) <- names
  colnames(coords) <- paste0("PC", seq_len(ncol(coords)) - 1)
  percent <- 100 * p$values$Relative_eig[positive]
  negative_percent <- -100 * sum(p$values$Relative_eig[negative])
  cat(name, ": ", n, " individuals, ", sum(positive), " positive eigenvalues, ",
      sum(negative), " negative, negative eigenvalues percent ", negative_percent, "\n", sep = "")
  write.table(format(coords, digits = 15), file.path(out, paste0(name, ".pcoa.r.projections.tsv")),
              sep = "\t", quote = FALSE)
  write.table(format(percent, digits = 15), file.path(out, paste0(name, ".pcoa.r.percent.tsv")),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(format(negative_percent, digits = 15), file.path(out, paste0(name, ".pcoa.r.negative.tsv")),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

gd_kosman <- function(name) {
  read.table(file.path("tests", "reference", "dists", paste0(name, ".gdkosman.tsv")), header = TRUE)$dist
}

write_pcoa("small", c(0.2, 0.3, 0.9, 0.9, 0.1, 0.8, 0.7, 0.7, 0.8, 0.2), paste0("i", 1:5))
write_pcoa("panel", gd_kosman("panel"), sprintf("s%03d", 0:199))
write_pcoa("four_alleles", gd_kosman("four_alleles"), sprintf("i%02d", 0:39))
