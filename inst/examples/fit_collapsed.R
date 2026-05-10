## Minimal example: fit the collapsed cgeneric on simulated genomic-style
## blocks. Mirrors fit_blocks.R, using the gmrfs package.

Sys.setenv(MKL_THREADING_LAYER = "SEQUENTIAL",
           MKL_NUM_THREADS     = "1",
           OMP_NUM_THREADS     = "1")

library(INLA)
library(Matrix)
library(gmrfs)

set.seed(42)
n <- 100; p <- 10; m_features <- 1000
freq_range <- c(0.05, 0.5)
sigma_u2   <- c(1.5, 1, 0.5, rep(0, p - 3))
sigma_e2   <- 0.3

G_list <- vector("list", p); Z_list <- G_list
for (j in seq_len(p)) {
  af  <- runif(m_features, freq_range[1], freq_range[2])
  raw <- matrix(rbinom(n * m_features, 2, rep(af, each = n)), n, m_features)
  Z   <- sweep(sweep(raw, 2, 2 * af, `-`), 2, sqrt(2 * af * (1 - af)), `/`)
  Gj  <- tcrossprod(Z) / m_features
  Gj  <- (Gj + t(Gj)) / 2 + 1e-8 * diag(n)
  Z_list[[j]] <- Z; G_list[[j]] <- Gj
}
u_list <- lapply(seq_len(p), function(j) {
  if (sigma_u2[j] == 0) numeric(n)
  else as.numeric(sqrt(sigma_u2[j] / m_features) *
                    (Z_list[[j]] %*% rnorm(m_features)))
})
y <- Reduce(`+`, u_list) + rnorm(n, sd = sqrt(sigma_e2))

cm  <- cblocks(n = n, p = p, G = do.call(rbind, G_list))
idx <- seq_len(n)

fit <- inla(y ~ 1 + f(idx, model = cm),
            data            = data.frame(y = y, idx = idx),
            control.family  = list(hyper = list(prec = list(prior = "pc.prec",
                                                            param = c(10, 0.01)))),
            control.compute = list(config = TRUE),
            control.inla    = list(int.strategy = "eb"),
            num.threads     = "1:1")
print(fit$summary.hyperpar)
