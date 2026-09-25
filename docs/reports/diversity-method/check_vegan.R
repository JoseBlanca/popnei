# The rarefied number of alleles of each population of the panel, from vegan.
suppressMessages(library(vegan))
cat("vegan", as.character(packageVersion("vegan")), "\n")
d <- read.table("panel_counts.tsv", header = TRUE)
g <- 20
for (p in sort(unique(d$pop))) {
  s <- d[d$pop == p, ]
  s <- s[s$called >= g, ]
  m <- as.matrix(s[, c("n0", "n1")])
  e <- rarefy(m, sample = g)
  cat(sprintf("%s  variants %d  mean expected alleles in a draw of %d: %.10f\n",
              p, nrow(s), g, mean(e)))
}
