S <- Sys.getenv("S"); .libPaths(c(file.path(S,"rlib"), .libPaths()))
suppressMessages(library(adegenet)); source(file.path(S,"pgr/PopGenReport/R/gd.kosman.r"))
for (name in c("tetra","haploid","tworked","hworked")) {
  df <- read.csv(file.path(S,paste0(name,".gts.csv")), row.names=1, colClasses="character", na.strings="NA")
  k <- if (startsWith(name,"t")) 4 else 1
  gi <- df2genind(df, sep="/", ploidy=k, NA.char="NA", type="codom")
  r <- gd.kosman(gi); d <- r$geneticdist; n <- r$loci_used; ni <- nrow(d); vec <- c(); nv <- c()
  for (i in 1:(ni-1)) for (j in (i+1):ni) { vec <- c(vec, d[j,i]); nv <- c(nv, n[j,i]) }
  if (file.exists(file.path(S,paste0(name,".python.txt")))) { py <- scan(file.path(S,paste0(name,".python.txt")), quiet=TRUE); cat(name, "pairs", length(vec), "max abs diff R vs python", max(abs(py-vec)), "\n") }
  else { cat(name, "\n"); print(vec); print(nv) }
  if (name %in% c("tetra","haploid")) cat("  first", format(vec[1:3], digits=17), " n", nv[1:3], "\n")
}
