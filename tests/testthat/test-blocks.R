ar1_cov <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
n <- 8
G_list <- list(ar1_cov(n, 0.5), ar1_cov(n, 0.2), ar1_cov(n, 0.8))
theta <- c(0.3, -0.5, 1.2)
Sigma <- function(th, G) Reduce(`+`, Map(function(t, g) exp(-t) * g, th, G))

test_that("cblocks has covariance sum_j exp(-theta_j) G_j", {
  cm <- cblocks(n, 3, G_list)
  S <- Sigma(theta, G_list)
  expect_equal(as.matrix(solve(gmrf_Q(cm, theta))), S, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(gmrf_log_norm_const(cm, theta),
               -0.5 * as.numeric(determinant(S)$modulus) - 0.5 * n * log(2 * pi))
  expect_equal(gmrf_initial(cm), rep(1, 3))
  expect_equal(gmrf_mu(cm, theta), numeric(n))
  expect_true(is.finite(gmrf_log_prior(cm, theta)))
  expect_equal(dim(gmrf_graph(cm)), c(n, n))
})

test_that("theta is checked before it reaches the C code", {
  cm <- cblocks(n, 3, G_list)
  expect_error(gmrf_Q(cm, c(0, 0)), "length 3")
  expect_error(gmrf_Q(cm, c(0, NA, 0)), "missing")
  expect_error(gmrf_Q(list()), "model handle")
})

test_that("fixed precisions leave no hyperparameters", {
  cm <- cblocks(n, 3, G_list, fixed = exp(theta))
  expect_length(gmrf_initial(cm), 0)
  expect_equal(as.matrix(solve(gmrf_Q(cm))), Sigma(theta, G_list),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("sblocks is block-diagonal with blocks exp(theta_j) R_j", {
  R_list <- lapply(G_list, solve)
  sm <- sblocks(n, 3, R_list)
  Q <- as.matrix(gmrf_Q(sm, theta))
  ref <- as.matrix(Matrix::bdiag(Map(function(t, r) exp(t) * r, theta, R_list)))
  expect_equal(Q, ref, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(gmrf_log_norm_const(sm, theta),
               0.5 * as.numeric(determinant(ref)$modulus) - 0.5 * 3 * n * log(2 * pi))
})

test_that("pblocks separates one block and collapses the rest", {
  pm <- pblocks(n, 3, G_list, separate = 2)
  expect_equal(attr(pm, "block_order"), c(2L, 1L, 3L))
  Q <- as.matrix(gmrf_Q(pm, theta))
  ## theta[1] is the separated block, theta[2:3] the others in order
  expect_equal(Q[1:n, 1:n], exp(theta[1]) * solve(G_list[[2]]),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(solve(Q[n + 1:n, n + 1:n]), Sigma(theta[2:3], G_list[c(1, 3)]),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(max(abs(Q[1:n, n + 1:n])), 0)
})

test_that("omp_select finds the active region", {
  set.seed(1)
  m <- 100
  Z <- lapply(1:8, function(j) matrix(rnorm(m * 5), m, 5))
  y <- Z[[3]] %*% rnorm(5) + rnorm(m, sd = 0.5)
  sel <- omp_select(y, lapply(Z, tcrossprod))
  expect_s3_class(sel, "omp_select")
  expect_equal(sel$selected[1], 3)
})
