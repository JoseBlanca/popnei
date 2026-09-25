# The curve of expected r2 fitted to the pairs decay.py wrote: the table of
# "How it is verified" of the item "LD against distance, per population" of
# docs/specs/ld.md, and the numbers its text quotes.
W <- Sys.getenv("LD_WORK", ".")

# Hill and Weir (1988) with the correction for the n gametes sampled, and
# the same curve with that correction dropped.
hw <- function(d, C, n) { p <- C * d
  ((10 + p) / ((2 + p) * (11 + p))) *
    (1 + ((3 + p) * (12 + 12 * p + p * p)) / (n * (2 + p) * (11 + p))) }
plain <- function(d, C, n) { p <- C * d; (10 + p) / ((2 + p) * (11 + p)) }
sved <- function(d, C, n) 1 / (1 + C * d)

# The sum popnei minimises: the pairs at a distance enter through their
# number and the mean of their r2 alone.
ss_of <- function(d, w, y, f, n) function(C) sum(w * (y - f(d, C, n))^2)
brent <- function(d, w, y, f, n, lo = -12, hi = 2) {
  o <- optimize(function(l) ss_of(d, w, y, f, n)(10^l), c(lo, hi), tol = 1e-14)
  list(C = 10^o$minimum, ss = o$objective) }
gauss_newton <- function(d, w, y, f, n, start) {
  t <- data.frame(d = d, y = y, w = w)
  unname(coef(nls(y ~ f(d, C, n), data = t, weights = t$w, start = list(C = start),
                  algorithm = "port", lower = 1e-14, upper = 1e2,
                  control = nls.control(tol = 1e-10, maxiter = 500)))["C"]) }
# The curve falls without turning, so one rho gives half of its value at 0.
rho_at_half <- function(f, n)
  uniroot(function(p) f(p, 1, n) - f(0, 1, n) / 2, c(0, 1e6), tol = 1e-14)$root

read_pairs <- function(tag) {
  t <- read.table(file.path(W, paste0(tag, ".decay.tsv")), header = TRUE)
  t$mean_r2 <- t$sum_r2 / t$num_pairs; t }

# n is the individuals of the population, which decay_truth.py measured
# to be what the sample term of the curve is about.
pops <- list(all = 100, pop_a = 50, pop_b = 50)
cat("the table of the spec: the curve fitted to every pair\n")
fits <- list()
for (tag in names(pops)) {
  n <- pops[[tag]]; t <- read_pairs(tag)
  b <- brent(t$dist, t$num_pairs, t$mean_r2, hw, n)
  g <- gauss_newton(t$dist, t$num_pairs, t$mean_r2, hw, n, b$C)
  fits[[tag]] <- b
  cat(sprintf("  %-6s n=%3d  %d distances, %d pairs\n", tag, n, nrow(t), sum(t$num_pairs)))
  cat(sprintf("    rho per bp %.17g   r2 at 0 %.17g   half dist %.17g\n",
              b$C, hw(0, b$C, n), rho_at_half(hw, n) / b$C))
  cat(sprintf("    nls gives %.17g, which is %.3g of itself away from optimize\n",
              g, abs(g / b$C - 1)))
}

t <- read_pairs("all"); n <- 100; b <- fits[["all"]]
cat("\nhow well the curve describes this dataset, at the ten bins of the spec\n")
k <- pmin(floor((t$dist - 1) / 25000), 9)
np <- tapply(t$num_pairs, k, sum); sr <- tapply(t$sum_r2, k, sum)
mid <- 1 + (as.numeric(names(np)) + 0.5) * 25000
for (i in seq_along(np))
  cat(sprintf("  middle of the bin %7.0f  mean r2 %.4f  the curve %.4f\n",
              mid[i], sr[i] / np[i], hw(mid[i], b$C, n)))

cat("\nSved's curve on the same pairs, which is not the one used\n")
s <- brent(t$dist, t$num_pairs, t$mean_r2, sved, n)
cat(sprintf("  Hill and Weir leave the sum at %.2f, Sved at %.2f\n", b$ss, s$ss))
cat(sprintf("  Sved: rho per bp %.10g, r2 at 0 %g, half dist %.10g\n",
            s$C, sved(0, s$C, n), rho_at_half(sved, n) / s$C))

cat("\nwhat fitting the bin means instead of the pairs would give\n")
half <- rho_at_half(hw, n) / b$C
for (nb in c(10, 50)) {
  w <- 250000 / nb
  k <- pmin(floor((t$dist - 1) / w), nb - 1)
  np <- as.numeric(tapply(t$num_pairs, k, sum)); sr <- as.numeric(tapply(t$sum_r2, k, sum))
  mid <- 1 + (as.numeric(names(tapply(t$num_pairs, k, sum))) + 0.5) * w
  c2 <- brent(mid, np, sr / np, hw, n)
  cat(sprintf("  %2d bins at their middle: half dist %.10g, %.3g of the pairs' %.10g away\n",
              nb, rho_at_half(hw, n) / c2$C, rho_at_half(hw, n) / c2$C / half - 1, half))
}

cat("\nwhat the other answers to n would have given for the first row\n")
for (nn in c(100, 200)) {
  c2 <- brent(t$dist, t$num_pairs, t$mean_r2, hw, nn)
  cat(sprintf("  n = %3d: half dist %.10g\n", nn, rho_at_half(hw, nn) / c2$C))
}
c3 <- brent(t$dist, t$num_pairs, t$mean_r2, plain, n)
cat(sprintf("  no sample term: half dist %.10g\n", rho_at_half(plain, n) / c3$C))

cat("\nwhat the ends of the searched range do to the first row\n")
for (lo in c(-12, -9)) for (hi in c(2, 0)) {
  r <- brent(t$dist, t$num_pairs, t$mean_r2, hw, n, lo, hi)
  cat(sprintf("  10^[%3d, %2d]: rho per bp %.17g, %.3g of itself away\n",
              lo, hi, r$C, abs(r$C / b$C - 1)))
}

cat("\nthe rho at which the curve is half of its value at 0, by n\n")
for (nn in c(25, 50, 100, 200))
  cat(sprintf("  n = %3d: r2 at 0 %.17g, rho at half %.17g\n",
              nn, hw(0, 1, nn), rho_at_half(hw, nn)))

cat("\nthe case of the first cargo test: r2 taken from the curve itself\n")
C0 <- 1e-4; d0 <- seq(1000, 250000, by = 1000); y0 <- hw(d0, C0, 100)
r0 <- brent(d0, rep(1, length(d0)), y0, hw, 100)
cat(sprintf("  made at rho per bp %g, recovered %.17g, %.3g of itself away\n",
            C0, r0$C, abs(r0$C / C0 - 1)))
cat(sprintf("  half dist %.17g\n", rho_at_half(hw, 100) / C0))
cat(sprintf("  r2 at 1000, 2000, 3000 bp: %.17g %.17g %.17g\n", y0[1], y0[2], y0[3]))
cat(sprintf("\nR %s.%s\n", R.version$major, R.version$minor))
