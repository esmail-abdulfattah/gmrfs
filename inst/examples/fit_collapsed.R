## Minimal smoke test for gmrfs.
##
## Run with:
##   Rscript -e 'source(system.file("examples/fit_collapsed.R", package = "gmrfs"))'
## or directly:
##   Rscript inst/examples/fit_collapsed.R
##
## Builds a small genomic-style block model, fits with `cblocks()`, and
## prints the hyperparameter summary.  Should finish in a few seconds.

suppressPackageStartupMessages({
  library(INLA); library(Matrix); library(gmrfs)
})

set.seed(42)
n <- 100; p <- 3; m_features <- 500
freq_range <- c(0.05, 0.5)
sigma_u2   <- c(1.5, 1.0, 0.5)
sigma_e2   <- 0.3

G_list <- vector("list", p); Z_list <- G_list
for (j in seq_len(p)) {
  af  <- runif(m_features, freq_range[1], freq_range[2])
  raw <- matrix(rbinom(n * m_features, 2, rep(af, each = n)), n, m_features)
  Z   <- sweep(sweep(raw, 2, 2 * af, `-`), 2, sqrt(2 * af * (1 - af)), `/`)
  Gj  <- tcrossprod(Z) / m_features
  Z_list[[j]] <- Z
  G_list[[j]] <- (Gj + t(Gj)) / 2 + 1e-8 * diag(n)
}
u_list <- lapply(seq_len(p), function(j)
  as.numeric(sqrt(sigma_u2[j] / m_features) *
             (Z_list[[j]] %*% rnorm(m_features))))
y <- Reduce(`+`, u_list) + rnorm(n, sd = sqrt(sigma_e2))

cm  <- cblocks(n = n, p = p, G_list = G_list, tau0 = 1.0)
fit <- inla(y ~ 1 + f(idx, model = cm),
            data            = data.frame(y = y, idx = seq_len(n)),
            control.family  = list(hyper = list(prec = list(prior = "pc.prec",
                                                            param = c(10, 0.01)))),
            control.compute = list(config = TRUE),
            control.inla    = list(int.strategy = "eb"))

cat("\n[gmrfs smoke test] OK\n\n")
print(fit$summary.hyperpar)
