# What ape's pcoa gives for five distance matrices, the literals of the
# principal coordinates of docs/specs/pca.md. Run from the root of the
# repository with R 4.6.1 and ape 5.8.1:
#
#     Rscript tests/reference/pca/pcoa_reference.R tests/reference/pca
#
# The five matrices:
#
# - small, the ten distances of five individuals of test_pcoa of pyNei's
#   test/test_pca.py, which are not Euclidean: one eigenvalue is negative.
# - small_twin, small with a sixth individual, i6, at distance 0 from i5 and
#   at the distances of i5 from the others, so that the eigenvalue 0 of the
#   centering has a second eigenvector beside it.
# - panel, the Kosman distances of the 200 individuals of
#   tests/reference/dists/panel.vcf.gz, as R's gd.kosman gave them in
#   tests/reference/dists/panel.gdkosman.tsv.
# - panel_clone, the panel with a 201st individual, s200, whose genotypes are
#   those of s000, as make_panel_clone.py writes it into panel_clone.vcf.gz:
#   its Kosman distances are those of the panel with s200 at 0 from s000 and
#   at the distances of s000 from the others. Its eigenvalue 0 has two
#   eigenvectors, as small_twin's has.
# - four_alleles, the Kosman distances of the 40 individuals of
#   tests/reference/dists/four_alleles.vcf.gz, from four_alleles.gdkosman.tsv,
#   which have no negative eigenvalue.
#
# An eigenvalue counts as positive or negative when it is further from 0
# than the largest times n times the machine epsilon, the threshold of the
# spec. Each component gets the sign of the spec: its projection of the
# largest absolute value is positive, and among projections within 64 units
# in the last place of that value, the first individual's.
#
# For each matrix <name> it writes, when the matrix is Euclidean, what pcoa
# gives with no correction:
#
# - <name>.pcoa.r.projections.tsv, one row per individual, one column per
#   component;
# - <name>.pcoa.r.percent.tsv, each eigenvalue over the sum of all of them,
#   times 100, Relative_eig;
#
# and, for every matrix, what pcoa gives with correction = "lingoes", which
# adds 2c to every squared distance but those of an individual with itself,
# c being the absolute value of the most negative eigenvalue, 0 for a
# Euclidean matrix:
#
# - <name>.lingoes.r.projections.tsv and <name>.lingoes.r.percent.tsv, the
#   same of the corrected matrix, from vectors.cor and Rel_corr_eig;
# - <name>.lingoes.r.constant.tsv, c, and on a second line the negative
#   eigenvalues percent of the matrix before the correction, 100 times the
#   sum of the absolute values of the negative eigenvalues over the sum of
#   all of them.
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

with_the_sign_rule <- function(coords) {
  for (k in seq_len(ncol(coords))) {
    largest <- max(abs(coords[, k]))
    first <- which(abs(coords[, k]) >= largest - 64 * .Machine$double.eps * largest)[1]
    coords[, k] <- coords[, k] * sign(coords[first, k])
  }
  coords
}

write_table <- function(x, name, what, names = FALSE) {
  write.table(format(x, digits = 15), file.path(out, paste0(name, what)),
              sep = "\t", quote = FALSE, col.names = names, row.names = names)
}

write_pcoa <- function(name, dist_vector, names) {
  d <- square_of(dist_vector)
  n <- nrow(d)
  plain <- pcoa(as.dist(d), correction = "none")
  eig <- plain$values$Eigenvalues
  tolerance <- max(eig) * n * .Machine$double.eps
  positive <- eig > tolerance
  negative <- eig < -tolerance
  negative_percent <- -100 * sum(plain$values$Relative_eig[negative])
  if (any(negative)) {
    constant <- abs(min(eig))
    corrected <- pcoa(as.dist(d), correction = "lingoes")
    corr_eig <- corrected$values$Corr_eig
    coords <- corrected$vectors.cor
    percent <- 100 * corrected$values$Rel_corr_eig[corr_eig > max(corr_eig) * n * .Machine$double.eps]
  } else {
    constant <- 0
    coords <- plain$vectors[, seq_len(sum(positive)), drop = FALSE]
    percent <- 100 * plain$values$Relative_eig[positive]
    kept <- with_the_sign_rule(coords)
    dimnames(kept) <- list(names, paste0("PC", seq_len(ncol(kept)) - 1))
    write_table(kept, name, ".pcoa.r.projections.tsv", TRUE)
    write_table(percent, name, ".pcoa.r.percent.tsv")
  }
  coords <- with_the_sign_rule(coords[, seq_along(percent), drop = FALSE])
  dimnames(coords) <- list(names, paste0("PC", seq_len(ncol(coords)) - 1))
  cat(name, ": ", n, " individuals, ", sum(positive), " positive eigenvalues, ",
      sum(negative), " negative, negative eigenvalues percent ", negative_percent,
      ", Lingoes constant ", constant, ", ", ncol(coords), " components after it\n", sep = "")
  write_table(coords, name, ".lingoes.r.projections.tsv", TRUE)
  write_table(percent, name, ".lingoes.r.percent.tsv")
  write_table(c(constant, negative_percent), name, ".lingoes.r.constant.tsv")
}

gd_kosman <- function(name) {
  read.table(file.path("tests", "reference", "dists", paste0(name, ".gdkosman.tsv")), header = TRUE)$dist
}

small <- c(0.2, 0.3, 0.9, 0.9, 0.1, 0.8, 0.7, 0.7, 0.8, 0.2)
write_pcoa("small", small, paste0("i", 1:5))
# i6 is i5 again: its distances to i1 .. i4 are those of i5, and 0 to i5. In
# the order of the vector the pairs of i6 come last of each row.
twin <- square_of(small)
twin <- rbind(cbind(twin, twin[, 5]), c(twin[5, ], 0))
write_pcoa("small_twin", twin[lower.tri(twin)], paste0("i", 1:6))
write_pcoa("panel", gd_kosman("panel"), sprintf("s%03d", 0:199))
panel <- square_of(gd_kosman("panel"))
panel_clone <- rbind(cbind(panel, panel[, 1]), c(panel[1, ], 0))
write_pcoa("panel_clone", panel_clone[lower.tri(panel_clone)], sprintf("s%03d", 0:200))
write_pcoa("four_alleles", gd_kosman("four_alleles"), sprintf("i%02d", 0:39))
