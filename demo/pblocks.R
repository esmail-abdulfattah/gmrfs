suppressPackageStartupMessages({
  library(gmrfs); library(INLA); library(Matrix)
})

set.seed(1)
n <- 20; p <- 2

ar1_cov <- function(n, rho) outer(seq_len(n), seq_len(n),
                                  function(i, j) rho^abs(i - j))
G_list  <- list(ar1_cov(n, 0.5), ar1_cov(n, 0.2))

## pblocks: latent dim 2n. Block 1 separate, block 2 collapsed-into-rest.
cm  <- pblocks(n = n, p = p, G_list = G_list, separate = 1L, tau0 = 1.0)

## Each obs i maps to latent positions (i, n+i).
y   <- rnorm(n, 0, 0.5)
A   <- sparseMatrix(i = rep(seq_len(n), 2),
                    j = c(seq_len(n), n + seq_len(n)),
                    x = 1.0, dims = c(n, 2L * n))
stk <- inla.stack(data    = list(y = y),
                  A       = list(A, 1),
                  effects = list(idx = seq_len(2L * n),
                                 intercept = rep(1, n)),
                  tag     = "obs")

m <- INLA::inla(y ~ -1 + intercept + f(idx, model = cm),
                data              = inla.stack.data(stk),
                control.predictor = list(A = inla.stack.A(stk)),
                family            = "gaussian",
                control.inla      = list(int.strategy = "eb"))

print(m$summary.hyperpar)
stopifnot(isTRUE(m$ok))
cat("\n\n[OK] pblocks minimal demo ran successfully.\n")
