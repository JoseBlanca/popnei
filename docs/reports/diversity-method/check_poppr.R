# The alleles each population of the panel called, and the private ones.
suppressMessages({library(adegenet); library(poppr)})
cat("adegenet", as.character(packageVersion("adegenet")),
    " poppr", as.character(packageVersion("poppr")), "\n")
d <- read.table("panel_genotypes.tsv", header = TRUE, sep = "\t",
                colClasses = "character", check.names = FALSE)
pop <- factor(d$pop)
gt <- d[, setdiff(names(d), c("ind", "pop"))]
rownames(gt) <- d$ind
ind <- df2genind(gt, sep = "/", ploidy = 2, pop = pop, NA.char = "NA")
tab <- tab(genind2genpop(ind, quiet = TRUE))
cat("\nalleles called by each population, over the 1200 variants\n")
print(rowSums(tab > 0))
pa <- private_alleles(ind, count.alleles = FALSE)
cat("\nprivate alleles of each population\n")
print(rowSums(pa))
cat("\nexpected heterozygosity per population, adegenet Hs (no sample correction)\n")
print(Hs(ind))
