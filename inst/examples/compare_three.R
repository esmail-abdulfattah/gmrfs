## compare_three.R -- recover sigma_u2 from one simulated y using all
## three block formulations: collapsed, partial, stacked.
##
##   y_i = sum_{j=1}^p u_j[i] + eps_i
##   u_j ~ N(0, sigma_u2[j] * G_j),   eps ~ N(0, sigma_e2)
##
## All three cgeneric models target the same posterior for the per-block
## log-precisions theta_j; they differ in the latent parameterisation
## (and therefore in cost / numerics):
##
##   collapsed   latent dim n     Q = (sum_j e^{-theta_j} G_j)^{-1}
##   partial     latent dim 2n    block 1 split, blocks 2..p collapsed
##   stacked     latent dim p*n   each block separate, block-diagonal Q

suppressPackageStartupMessages({
  library(INLA); library(Matrix); library(gmrfs)
})

## ============================================================
## 1. Simulate (small enough to run in seconds)
## ============================================================
set.seed(42)
n          <- 100
p          <- 3
m_features <- 500
freq_range <- c(0.05, 0.5)
sigma_u2   <- c(1.5, 1.0, 0.5)
sigma_e2   <- 0.3
stopifnot(length(sigma_u2) == p)

G_list <- Z_list <- vector("list", p)
for (j in seq_len(p)) {
  af  <- runif(m_features, freq_range[1], freq_range[2])
  raw <- matrix(rbinom(n * m_features, 2, rep(af, each = n)), n, m_features)
  Z   <- sweep(sweep(raw, 2, 2 * af, `-`), 2, sqrt(2 * af * (1 - af)), `/`)
  Gj  <- tcrossprod(Z) / m_features
  Z_list[[j]] <- Z
  G_list[[j]] <- (Gj + t(Gj)) / 2 + 1e-8 * diag(n)
}
u_list <- lapply(seq_len(p), function(j) {
  if (sigma_u2[j] == 0) numeric(n)
  else as.numeric(sqrt(sigma_u2[j] / m_features) *
                  (Z_list[[j]] %*% rnorm(m_features)))
})
y <- Reduce(`+`, u_list) + rnorm(n, sd = sqrt(sigma_e2))

## per-block precision matrices (needed by stacked)
R_list <- lapply(G_list, function(Gj) chol2inv(chol(Gj)))

## ============================================================
## 2. Common INLA settings
## ============================================================
inla_args <- list(
  control.family  = list(hyper = list(prec = list(prior = "pc.prec",
                                                  param = c(10, 0.01)))),
  control.compute = list(config = TRUE),
  control.inla    = list(int.strategy = "eb"),
  verbose         = FALSE
)

## Posterior mean of sigma_j^2 = exp(-theta_j) via the marginal of theta_j.
sigma2_mean <- function(fit, j) {
  m <- inla.tmarginal(function(t) exp(-t),
                      fit$marginals.hyperpar[[j + 1L]])
  inla.emarginal(function(x) x, m)
}

summarize <- function(name, fit) {
  cat(sprintf("\n=========  %s  =========\n", name))
  est <- vapply(seq_len(p), function(j) sigma2_mean(fit, j), numeric(1))
  out <- data.frame(
    region      = paste0("region_", seq_len(p)),
    true_sigma2 = round(sigma_u2, 3),
    post_mean   = round(est,      3)
  )
  print(out)
  cat(sprintf("CPU = %.1fs\n", as.numeric(fit$cpu.used[["Total"]])))
  est
}

## ============================================================
## 3. Collapsed:  latent dim = n,  f(idx) directly
## ============================================================
cm_col  <- cblocks(n = n, p = p, G_list = G_list, tau0 = 1.0)
fit_col <- do.call(inla, c(
  list(formula = y ~ 1 + f(idx, model = cm_col),
       data    = data.frame(y = y, idx = seq_len(n))),
  inla_args))
est_col <- summarize("collapsed", fit_col)

## ============================================================
## 4. Partial:  latent dim = 2n,  A maps obs i -> (i, n+i)
## ============================================================
cm_par  <- pblocks(n = n, p = p, G_list = G_list, separate = 1L, tau0 = 1.0)
N_par   <- 2L * n
A_par   <- sparseMatrix(i    = rep(seq_len(n), 2),
                        j    = c(seq_len(n), n + seq_len(n)),
                        x    = 1.0,
                        dims = c(n, N_par))
stk_par <- inla.stack(data    = list(y = y),
                      A       = list(A_par, 1),
                      effects = list(idx       = seq_len(N_par),
                                     intercept = rep(1, n)),
                      tag     = "obs")
fit_par <- do.call(inla, c(
  list(formula = y ~ -1 + intercept + f(idx, model = cm_par),
       data    = inla.stack.data(stk_par),
       control.predictor = list(A = inla.stack.A(stk_par))),
  inla_args))
est_par <- summarize("partial", fit_par)

## ============================================================
## 5. Stacked:  latent dim = p*n, A maps obs i -> (i, n+i, ..., (p-1)*n+i)
## ============================================================
cm_stk  <- sblocks(n = n, p = p, R_list = R_list, tau0 = 1.0)
N_stk   <- p * n
A_stk   <- sparseMatrix(i    = rep(seq_len(n), p),
                        j    = as.vector(outer(seq_len(n), (0:(p - 1)) * n, `+`)),
                        x    = 1.0,
                        dims = c(n, N_stk))
stk_stk <- inla.stack(data    = list(y = y),
                      A       = list(A_stk, 1),
                      effects = list(idx       = seq_len(N_stk),
                                     intercept = rep(1, n)),
                      tag     = "obs")
fit_stk <- do.call(inla, c(
  list(formula = y ~ -1 + intercept + f(idx, model = cm_stk),
       data    = inla.stack.data(stk_stk),
       control.predictor = list(A = inla.stack.A(stk_stk))),
  inla_args))
est_stk <- summarize("stacked", fit_stk)

## ============================================================
## 6. Side-by-side: posterior mean of sigma_j^2
## ============================================================
cmp <- data.frame(
  region      = paste0("region_", seq_len(p)),
  true_sigma2 = round(sigma_u2, 3),
  collapsed   = round(est_col,  3),
  partial     = round(est_par,  3),
  stacked     = round(est_stk,  3)
)

cat("\n=========  Posterior mean of sigma_j^2  =========\n")
print(cmp)
cat(sprintf("\nCPU times  collapsed=%.1fs   partial=%.1fs   stacked=%.1fs\n",
            as.numeric(fit_col$cpu.used[["Total"]]),
            as.numeric(fit_par$cpu.used[["Total"]]),
            as.numeric(fit_stk$cpu.used[["Total"]])))
