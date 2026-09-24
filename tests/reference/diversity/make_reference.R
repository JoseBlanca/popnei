#!/usr/bin/env Rscript
#
# The alleles the populations of the panel call, from vegan, adegenet and poppr.
#
# It writes the reference numbers of three items of docs/specs/diversity.md,
# the number of alleles, the private alleles and the variable variants, which
# the tests of the diversity module compare against. Run it from anywhere:
#
#     Rscript tests/reference/diversity/make_reference.R
#
# It reads the panel of docs/specs/stats.md, tests/reference/stats/panel.vcf.gz
# with the populations of its individuals in
# tests/reference/stats/panel_pops_bcftools.txt: 1200 biallelic diploid
# variants of 200 individuals, 3 in 100 genotypes missing whole, in the
# populations p0, p1 and p2 of 48, 68 and 84 individuals. It finds both from
# its own location, and it reads no file another script writes: the genotype
# table adegenet takes and the per variant allele counts vegan takes are both
# built here out of that VCF.
#
# It writes three files beside itself, one per item of the spec, and the name
# of a column says which program gave the number in it:
#
#     panel_num_alleles.tsv      the alleles each population called, from
#                                adegenet, and the alleles a draw of 20 of its
#                                called alleles is expected to show, from vegan
#     panel_private_alleles.tsv  the alleles a population called that no other
#                                population called, from poppr
#     panel_variable_vars.tsv    the variants that vary in each population,
#                                from adegenet, and the chance that such a
#                                draw of 20 varies, from vegan
#
# The tables of the spec are computed with min_num_individuals 20, how many
# called genotypes a population needs at a variant for the variant to count
# for it. None of the three programs has such a threshold, and none is needed
# here: the fewest called genotypes any of the three populations has at any
# variant of the panel is 42, so all 1200 variants count for all three. The
# script stops if that stops holding, since then the programs would be
# counting variants that popnei leaves out.
#
# The standardized ratio of variable variants is what vegan's rarefy gives
# minus 1. On a variant of two alleles a draw shows one allele or two, so the
# alleles it is expected to show are one plus the chance that it varies. The
# script stops if any variant of the panel has more than two alleles, which is
# what that identity needs; "How it is verified" of the variable variants of
# the spec has the rest of it.
#
# The floats go in to 17 significant digits and not to the ten the spec
# prints, because the tests compare with vegan within 1e-12 of the value and a
# number rounded to ten decimals is up to 5e-11 away from it.
#
# It compares every number it got with the literal of the spec's tables and
# stops at the first one that differs, as
# tests/reference/stats/make_reference.py does for plink2 and bcftools, so the
# files it leaves hold the numbers of the spec or the script fails.
#
# The versions that wrote the committed files, run on 24 September 2026: R
# 4.6.1, vegan 2.7.6, adegenet 2.1.11 and poppr 2.9.8. The script stops on any
# other version of the three packages. hierfstat, whose allelic.richness
# computes the same rarefaction as vegan's rarefy, is not among them: it does
# not install on the owner's machine, RcppParallel and then its dependency
# gaston failing to build on 24 September 2026.

suppressMessages({
  library(vegan)
  library(adegenet)
  library(poppr)
})

VERSIONS <- c(vegan = "2.7.6", adegenet = "2.1.11", poppr = "2.9.8")

# num_called_alleles of the spec, the g every population is brought down to.
NUM_CALLED_ALLELES <- 20L

# min_num_individuals of the spec, in called genotypes.
MIN_NUM_INDIVIDUALS <- 20L

# The spec prints ten decimals, so its literal is up to 5e-11 from the value.
TOLERANCE <- 1e-10

check_versions <- function() {
  for (package in names(VERSIONS)) {
    found <- as.character(packageVersion(package))
    if (found != VERSIONS[[package]]) {
      stop(sprintf("%s %s is needed; found %s", package, VERSIONS[[package]], found))
    }
  }
}

script_dir <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  path <- sub("^--file=", "", arguments[startsWith(arguments, "--file=")])
  if (length(path) != 1) {
    stop("the script has to be run with Rscript, which passes --file=")
  }
  normalizePath(dirname(path))
}

# The genotypes of the panel as a matrix of variants by individuals, holding
# the two allele numbers of each call, with NA for a missing one.
read_panel <- function(vcf_path) {
  handle <- gzfile(vcf_path, "rt")
  lines <- readLines(handle)
  close(handle)
  header <- lines[startsWith(lines, "#CHROM")]
  if (length(header) != 1) {
    stop(sprintf("%s has no #CHROM header line", vcf_path))
  }
  individuals <- strsplit(header, "\t", fixed = TRUE)[[1]][-(1:9)]
  fields <- do.call(rbind, strsplit(lines[!startsWith(lines, "#")], "\t", fixed = TRUE))
  if (!all(fields[, 9] == "GT")) {
    stop(sprintf(paste("%s has a variant whose FORMAT is not GT alone, so a",
                       "call of it is not a pair of alleles"), vcf_path))
  }
  calls <- fields[, -(1:9), drop = FALSE]
  dimnames(calls) <- list(fields[, 3], individuals)
  if (!all(grepl("^([0-9]+|[.])[/|]([0-9]+|[.])$", calls))) {
    bad <- unique(calls[!grepl("^([0-9]+|[.])[/|]([0-9]+|[.])$", calls)])
    stop(sprintf("%s has calls that are not a pair of alleles: %s", vcf_path,
                 paste(utils::head(bad, 5), collapse = " ")))
  }
  list(
    individuals = individuals,
    calls = calls,
    first = allele_numbers(calls, 1L, vcf_path),
    second = allele_numbers(calls, 2L, vcf_path)
  )
}

# One side of every genotype as an allele number, NA where it is not called.
allele_numbers <- function(calls, side, vcf_path) {
  text <- if (side == 1L) sub("[/|].*$", "", calls) else sub("^.*[/|]", "", calls)
  numbers <- rep(NA_integer_, length(text))
  called <- text != "."
  numbers[called] <- suppressWarnings(as.integer(text[called]))
  if (anyNA(numbers[called])) {
    bad <- unique(text[called][is.na(numbers[called])])
    stop(sprintf("%s has allele numbers this script cannot read: %s",
                 vcf_path, paste(bad, collapse = " ")))
  }
  dim(numbers) <- dim(calls)
  numbers
}

read_pops <- function(pops_path, individuals) {
  table <- read.table(pops_path, sep = "\t", col.names = c("individual", "pop"),
                      colClasses = "character")
  if (!identical(table$individual, individuals)) {
    stop(sprintf("%s does not name the individuals of the VCF in their order", pops_path))
  }
  table$pop
}

# For each population, the times it called each allele at each variant, the
# called alleles it has there and the genotypes it called there.
count_alleles <- function(panel, pop_of, pop_names) {
  max_allele <- max(c(panel$first, panel$second), na.rm = TRUE)
  num_vars <- nrow(panel$calls)
  counts <- array(0L, c(num_vars, length(pop_names), max_allele + 1L))
  called_alleles <- matrix(0, num_vars, length(pop_names),
                           dimnames = list(NULL, pop_names))
  called_genotypes <- called_alleles
  for (index in seq_along(pop_names)) {
    columns <- which(pop_of == pop_names[index])
    first <- panel$first[, columns, drop = FALSE]
    second <- panel$second[, columns, drop = FALSE]
    for (allele in 0:max_allele) {
      counts[, index, allele + 1L] <- rowSums(first == allele, na.rm = TRUE) +
        rowSums(second == allele, na.rm = TRUE)
    }
    called_alleles[, index] <- rowSums(counts[, index, , drop = FALSE], dims = 2)
    called_genotypes[, index] <- rowSums(!is.na(first)) / 2 + rowSums(!is.na(second)) / 2
  }
  list(counts = counts, called_alleles = called_alleles,
       called_genotypes = called_genotypes, max_allele = max_allele)
}

# The four things the numbers below rest on: no genotype is half called, so
# adegenet and the counts built here read the same alleles; no variant has more
# than two alleles, which the standardized ratio of variable variants needs;
# and every variant reaches both thresholds in every population, so the
# programs, which have no thresholds, count the variants popnei counts.
check_the_panel <- function(panel, alleles, pop_names) {
  half_called <- xor(is.na(panel$first), is.na(panel$second))
  if (any(half_called)) {
    stop(sprintf(paste("half called genotypes in the panel: %d. adegenet and",
                       "the counts built here would read them differently"),
                 sum(half_called)))
  }
  if (alleles$max_allele != 1L) {
    stop(sprintf(paste("the panel has a variant of %d alleles; the standardized",
                       "ratio of variable variants is taken from vegan's rarefy,",
                       "which needs two"),
                 alleles$max_allele + 1L))
  }
  fewest <- apply(alleles$called_genotypes, 2, min)
  if (any(fewest < MIN_NUM_INDIVIDUALS)) {
    short <- pop_names[fewest < MIN_NUM_INDIVIDUALS]
    stop(sprintf(paste("min_num_individuals %d drops variants of %s, which the",
                       "programs here would go on counting"),
                 MIN_NUM_INDIVIDUALS, paste(short, collapse = " ")))
  }
  fewest_in_draw <- apply(alleles$called_alleles, 2, min)
  if (any(fewest_in_draw < NUM_CALLED_ALLELES)) {
    short <- pop_names[fewest_in_draw < NUM_CALLED_ALLELES]
    stop(sprintf(paste("a draw of %d called alleles leaves out variants of %s,",
                       "and the means written here are over every variant"),
                 NUM_CALLED_ALLELES, paste(short, collapse = " ")))
  }
}

# adegenet: the alleles each population called at each variant, as a table of
# one column per variant and allele holding how often the population called it.
run_adegenet <- function(panel, pop_of) {
  genotypes <- t(panel$calls)
  genotypes[startsWith(genotypes, ".")] <- "NA"
  individuals <- df2genind(as.data.frame(genotypes, stringsAsFactors = FALSE),
                           sep = "/", ploidy = 2, pop = factor(pop_of), NA.char = "NA")
  populations <- genind2genpop(individuals, quiet = TRUE)
  counts <- tab(populations)
  present <- counts > 0
  locus_of_column <- as.character(locFac(populations))
  per_variant <- t(apply(present, 1, function(row) tapply(row, locus_of_column, sum)))
  list(
    genind = individuals,
    alleles_called = rowSums(present),
    variable_vars = rowSums(per_variant > 1)
  )
}

# poppr: a table with a 1 where an allele is private to the population of the
# row, so the row sums are the private alleles of each population.
run_poppr <- function(genind, pop_names) {
  private <- private_alleles(genind, count.alleles = FALSE)
  totals <- setNames(rep(0L, length(pop_names)), pop_names)
  if (is.matrix(private) && ncol(private) > 0) {
    found <- rowSums(private)
    unknown <- setdiff(names(found), pop_names)
    if (length(unknown) > 0) {
      stop(sprintf("poppr named populations the panel does not have: %s",
                   paste(unknown, collapse = " ")))
    }
    totals[names(found)] <- as.integer(found)
  }
  totals
}

# vegan: the alleles a draw of NUM_CALLED_ALLELES of the called alleles of the
# population is expected to show, averaged over the variants.
run_vegan <- function(alleles, pop_names) {
  means <- setNames(rep(NA_real_, length(pop_names)), pop_names)
  kept <- setNames(rep(NA_integer_, length(pop_names)), pop_names)
  for (index in seq_along(pop_names)) {
    in_draw <- alleles$called_alleles[, index] >= NUM_CALLED_ALLELES
    counts <- matrix(alleles$counts[in_draw, index, ], nrow = sum(in_draw))
    expected <- rarefy(counts, sample = NUM_CALLED_ALLELES)
    means[index] <- mean(expected)
    kept[index] <- sum(in_draw)
  }
  list(alleles_in_draw = means, num_vars_in_draw = kept)
}

check_counts <- function(got, want, what) {
  if (!identical(as.numeric(unname(got)), as.numeric(want))) {
    stop(sprintf("%s: the spec has %s and the programs gave %s", what,
                 paste(want, collapse = " "), paste(unname(got), collapse = " ")))
  }
}

check_values <- function(got, want, what) {
  differences <- abs(unname(got) - want)
  if (any(differences > TOLERANCE)) {
    stop(sprintf("%s: the spec has %s and the programs gave %s, off by %s", what,
                 paste(sprintf("%.10f", want), collapse = " "),
                 paste(sprintf("%.10f", unname(got)), collapse = " "),
                 paste(sprintf("%.3g", differences), collapse = " ")))
  }
}

# The literals of the tables of docs/specs/diversity.md, for p0, p1 and p2.
check <- function(rows) {
  check_counts(rows$num_vars_with_data, c(1200L, 1200L, 1200L), "the variants that count")
  check_counts(rows$num_vars_in_draw, c(1200L, 1200L, 1200L), "the variants in the draw")
  check_counts(rows$alleles_called, c(2373L, 2377L, 2384L), "the alleles called")
  check_values(rows$alleles_mean, c(1.9775, 1.9808333333, 1.9866666667),
               "the mean alleles called")
  check_values(rows$alleles_in_draw, c(1.9283948650, 1.9219209943, 1.9197370844),
               "the alleles in a draw of 20")
  check_counts(rows$private_alleles, c(0L, 0L, 1L), "the private alleles")
  check_values(rows$private_mean, c(0, 0, 0.0008333333), "the mean private alleles")
  check_counts(rows$variable_vars, c(1173L, 1177L, 1184L), "the variable variants")
  check_values(rows$variable_ratio, c(0.9775, 0.9808333333, 0.9866666667),
               "the ratio of variable variants")
  check_values(rows$variable_in_draw, c(0.9283948650, 0.9219209943, 0.9197370844),
               "the ratio of variable variants in a draw of 20")
}

write_tsv <- function(path, header, columns) {
  lines <- header
  for (row in seq_along(columns[[1]])) {
    lines <- c(lines, paste(vapply(columns, function(column) as_field(column[[row]]), ""),
                            collapse = "\t"))
  }
  writeLines(lines, path)
}

as_field <- function(value) {
  if (is.character(value)) value
  else if (is.integer(value)) sprintf("%d", value)
  else sprintf("%.17g", value)
}

main <- function() {
  check_versions()
  here <- script_dir()
  stats <- file.path(dirname(here), "stats")
  panel <- read_panel(file.path(stats, "panel.vcf.gz"))
  pop_of <- read_pops(file.path(stats, "panel_pops_bcftools.txt"), panel$individuals)
  pop_names <- sort(unique(pop_of))
  alleles <- count_alleles(panel, pop_of, pop_names)
  check_the_panel(panel, alleles, pop_names)

  from_adegenet <- run_adegenet(panel, pop_of)
  from_poppr <- run_poppr(from_adegenet$genind, pop_names)
  from_vegan <- run_vegan(alleles, pop_names)
  if (!identical(names(from_adegenet$alleles_called), pop_names)) {
    stop(sprintf("adegenet gave the populations %s and the panel has %s",
                 paste(names(from_adegenet$alleles_called), collapse = " "),
                 paste(pop_names, collapse = " ")))
  }

  num_vars_with_data <- rep(nrow(panel$calls), length(pop_names))
  rows <- list(
    pop = pop_names,
    num_vars_with_data = as.integer(num_vars_with_data),
    num_vars_in_draw = as.integer(from_vegan$num_vars_in_draw),
    alleles_called = as.integer(from_adegenet$alleles_called),
    alleles_mean = from_adegenet$alleles_called / num_vars_with_data,
    alleles_in_draw = unname(from_vegan$alleles_in_draw),
    private_alleles = as.integer(from_poppr),
    private_mean = from_poppr / num_vars_with_data,
    variable_vars = as.integer(from_adegenet$variable_vars),
    variable_ratio = from_adegenet$variable_vars / num_vars_with_data,
    # One allele expected in a draw is one allele for certain plus the chance
    # of a second, which on a variant of two alleles is the chance it varies.
    variable_in_draw = unname(from_vegan$alleles_in_draw) - 1
  )
  check(rows)

  write_tsv(
    file.path(here, "panel_num_alleles.tsv"),
    "pop\tnum_vars_with_data\ttotal_adegenet\tmean\tnum_vars_in_draw\tin_draw_vegan",
    rows[c("pop", "num_vars_with_data", "alleles_called", "alleles_mean",
           "num_vars_in_draw", "alleles_in_draw")]
  )
  write_tsv(
    file.path(here, "panel_private_alleles.tsv"),
    "pop\tnum_vars_every_pop\ttotal_poppr\tmean",
    rows[c("pop", "num_vars_with_data", "private_alleles", "private_mean")]
  )
  write_tsv(
    file.path(here, "panel_variable_vars.tsv"),
    "pop\tnum_vars_with_data\ttotal_adegenet\tratio\tnum_vars_in_draw\tin_draw_vegan",
    rows[c("pop", "num_vars_with_data", "variable_vars", "variable_ratio",
           "num_vars_in_draw", "variable_in_draw")]
  )
  cat("done\n", file = stderr())
}

main()
