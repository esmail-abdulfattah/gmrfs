## compare_three.R -- recover sigma_u2 from one simulated y using all
## three formulations: collapsed, partial, stacked.
##
##   y_i = sum_{j=1}^p u_j[i] + eps_i
##   u_j ~ N(0, sigma_u2[j] * G_j),   eps ~ N(0, sigma_e2)
##
## All three cgeneric models target the same posterior for the per-
## block log-precisions theta_j = -log(sigma_u2[j]); they differ in
## the latent parameterization (and therefore in cost / numerics):
##
##   collapsed   latent dim n     Q = (sum_j e^{-theta_j} G_j)^{-1}
##   partial     latent dim 2n    block 1 split, blocks 2..p collapsed
##   stacked     latent dim p*n   each block separate, block-diagonal Q

Sys.setenv(MKL_THREADING_LAYER = "SEQUENTIAL",
           MKL_NUM_THREADS     = "1",
           OMP_NUM_THREADS     = "1")

suppressPackageStartupMessages({
  library(INLA); library(Matrix); library(gmrfs)
})

## ============================================================
## 1. Simulate
## ============================================================
set.seed(42)
n          <- 100
p          <- 5
m_features <- 1000
freq_range <- c(0.05, 0.5)
sigma_u2   <- c(1.5, 1.0, 0.5, 0.0, 0.0)   # last two regions are inactive
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

## per-block precision matrices (needed by partial and stacked)
R_list <- lapply(G_list, function(Gj) chol2inv(chol(Gj)))

true_theta <- ifelse(sigma_u2 > 0, -log(sigma_u2), NA_real_)

## ============================================================
## 2. Common INLA settings
## ============================================================
inla_args <- list(
  control.family  = list(hyper = list(prec = list(prior = "pc.prec",
                                                  param = c(10, 0.01)))),
  control.compute = list(config = TRUE),
  control.inla    = list(int.strategy = "eb"),
  num.threads     = "1:1",
  verbose         = FALSE
)

summarize <- function(name, fit) {
  hp        <- fit$summary.hyperpar
  theta_idx <- seq_len(p) + 1L                    # skip family precision
  out <- hp[theta_idx, c("mean", "sd",
                         "0.025quant", "0.5quant", "0.975quant", "mode")]
  rownames(out) <- paste0("theta_", seq_len(p))
  out$true  <- true_theta
  out$in_CI <- ifelse(is.na(true_theta), NA,
                      true_theta >= out[["0.025quant"]] &
                      true_theta <= out[["0.975quant"]])
  cat(sprintf("\n=========  %s  =========\n", name))
  print(round(out, 3))
  cat(sprintf("tau_y mode = %.3f  (truth = %.3f)\n",
              hp[1, "mode"], 1 / sigma_e2))
  cat(sprintf("CPU time   = %.1fs\n", as.numeric(fit$cpu.used[["Total"]])))
  out
}

## ============================================================
## 3. Collapsed:  latent dim = n,  f(idx) directly
## ============================================================
cm_col  <- cblocks(n = n, p = p, G = do.call(rbind, G_list))
fit_col <- do.call(inla, c(
  list(formula = y ~ 1 + f(idx, model = cm_col),
       data    = data.frame(y = y, idx = seq_len(n))),
  inla_args))
sum_col <- summarize("collapsed", fit_col)

## ============================================================
## 4. Partial:  latent dim = 2n,  A maps obs i -> (i, n+i)
##              x = (u_1, s_rest),  y_i = u_1[i] + s_rest[i] + eps_i
## ============================================================
cm_par <- pblocks(n = n, p = p,
                        G  = do.call(rbind, G_list),
                        R1 = R_list[[1]])
N_par <- 2L * n
A_par <- sparseMatrix(
  i    = rep(seq_len(n), 2),
  j    = c(seq_len(n), n + seq_len(n)),
  x    = 1.0,
  dims = c(n, N_par)
)
stk_par <- inla.stack(
  data    = list(y = y),
  A       = list(A_par, 1),
  effects = list(idx = seq_len(N_par), intercept = rep(1, n)),
  tag     = "obs"
)
fit_par <- do.call(inla, c(
  list(formula = y ~ -1 + intercept + f(idx, model = cm_par),
       data    = inla.stack.data(stk_par),
       control.predictor = list(A = inla.stack.A(stk_par))),
  inla_args))
sum_par <- summarize("partial", fit_par)

## ============================================================
## 5. Stacked:  latent dim = p*n, A maps obs i -> (i, n+i, ..., (p-1)*n+i)
##              x = (u_1, ..., u_p),  y_i = sum_j u_j[i] + eps_i
## ============================================================
R_blkdiag <- bdiag(R_list)                         # (p*n) x (p*n) sparse
cm_stk    <- sblocks(n = n, p = p, R = R_blkdiag)
N_stk     <- p * n
A_stk     <- sparseMatrix(
  i    = rep(seq_len(n), p),
  j    = as.vector(outer(seq_len(n), (0:(p - 1)) * n, `+`)),
  x    = 1.0,
  dims = c(n, N_stk)
)
stk_stk <- inla.stack(
  data    = list(y = y),
  A       = list(A_stk, 1),
  effects = list(idx = seq_len(N_stk), intercept = rep(1, n)),
  tag     = "obs"
)
fit_stk <- do.call(inla, c(
  list(formula = y ~ -1 + intercept + f(idx, model = cm_stk),
       data    = inla.stack.data(stk_stk),
       control.predictor = list(A = inla.stack.A(stk_stk))),
  inla_args))
sum_stk <- summarize("stacked", fit_stk)

## ============================================================
## 6. Side-by-side
## ============================================================
cmp <- data.frame(
  region    = paste0("theta_", seq_len(p)),
  true      = round(true_theta,    3),
  collapsed = round(sum_col$mode,  3),
  partial   = round(sum_par$mode,  3),
  stacked   = round(sum_stk$mode,  3)
)
cmp_sigma <- data.frame(
  region        = paste0("region_", seq_len(p)),
  true_sigma2   = round(sigma_u2, 3),
  collapsed     = round(exp(-sum_col$mode), 3),
  partial       = round(exp(-sum_par$mode), 3),
  stacked       = round(exp(-sum_stk$mode), 3)
)

cat("\n=========  Posterior mode of theta_j  =========\n")
print(cmp)
cat("\n=========  Implied sigma_j^2 = exp(-mode)  =========\n")
print(cmp_sigma)
cat(sprintf("\nCPU times  collapsed=%.1fs   partial=%.1fs   stacked=%.1fs\n",
            as.numeric(fit_col$cpu.used[["Total"]]),
            as.numeric(fit_par$cpu.used[["Total"]]),
            as.numeric(fit_stk$cpu.used[["Total"]])))
