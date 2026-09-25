.libPaths(c(file.path(Sys.getenv("S"),"rlib"), .libPaths()))
suppressMessages({library(adegenet); library(mmod); library(poppr)})
source(file.path(Sys.getenv("S"),"pgr/PopGenReport/R/gd.kosman.r"))
S <- Sys.getenv("S")
for (name in c("worked","sim_missing","multi")) {
  df <- read.csv(file.path(S,paste0(name,".gts.csv")), row.names=1, colClasses="character", na.strings="NA")
  gi <- df2genind(df, sep="/", ploidy=2, NA.char="NA", type="codom")
  k <- gd.kosman(gi)
  d <- k$geneticdist; n <- k$loci_used
  # condensed, row by row of the upper triangle: pairs (i,j) i<j ; the function fills the lower one
  ni <- nrow(d); vec <- c(); nvec <- c()
  for (i in 1:(ni-1)) for (j in (i+1):ni) { vec <- c(vec, d[j,i]); nvec <- c(nvec, n[j,i]) }
  write(format(vec, digits=17), file.path(S,paste0(name,".gdkosman.txt")), ncolumns=1)
  cat(name, "pairs", length(vec), "\n")
  if (name=="worked") { print(d); print(n); print(as.matrix(diss.dist(gi, percent=TRUE))); print(tryCatch(as.matrix(dist.codom(gi)), error=function(e) conditionMessage(e))) }
  if (name!="worked") { pn <- scan(file.path(S,paste0(name,".pynei.txt")), quiet=TRUE); cat(" max abs diff gd.kosman vs pyNei", max(abs(pn-vec)), "\n")
     dd <- as.matrix(diss.dist(gi, percent=TRUE)); v2 <- c(); for (i in 1:(ni-1)) for (j in (i+1):ni) v2 <- c(v2, dd[j,i]); cat(" max abs diff poppr diss.dist vs pyNei", max(abs(pn-v2)), "\n")
     cat(" first pairs", format(vec[1:3], digits=10), " n ", nvec[1:3], "\n") }
}
cat(R.version.string, as.character(packageVersion("adegenet")), as.character(packageVersion("poppr")), as.character(packageVersion("mmod")), "\n")
