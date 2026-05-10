suppressPackageStartupMessages({
  library(gmrfs); library(INLA); library(Matrix)
})

set.seed(1)
n <- 20; p <- 2

## Two structurally distinct block covariances (AR(1) with different rho).
ar1_cov <- function(n, rho) outer(seq_len(n), seq_len(n),
                                  function(i, j) rho^abs(i - j))
G_list  <- list(ar1_cov(n, 0.5), ar1_cov(n, 0.2))

cm  <- cblocks(n = n, p = p, G_list = G_list, tau0 = 1.0)

m <- INLA::inla(y ~ 1 + f(idx, model = cm),
                data   = list(y = rnorm(n, 0, 0.5), idx = seq_len(n)),
                family = "gaussian",
                control.inla = list(int.strategy = "eb"))

print(m$summary.hyperpar)
stopifnot(isTRUE(m$ok))
cat("\n\n[OK] cblocks minimal demo ran successfully.\n")
