ar1_cov <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
n <- 8
G_list <- list(ar1_cov(n, 0.5), ar1_cov(n, 0.2))
theta <- c(0.3, -0.5)

## Gaussian on A x = e, in orthonormal coordinates of that subspace
ld_sub <- function(Q, mu, A, e, x) {
  Q <- as.matrix(Q); q <- qr(t(A))
  N <- qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
  x0 <- as.numeric(t(A) %*% solve(tcrossprod(A), e))
  H <- t(N) %*% Q %*% N
  m <- as.numeric(-solve(H, t(N) %*% Q %*% (x0 - mu)))
  y <- as.numeric(t(N) %*% (x - x0)) - m
  -0.5 * ncol(N) * log(2 * pi) + 0.5 * as.numeric(determinant(H)$modulus) - 0.5 * sum(y * (H %*% y))
}

test_that("log density without constraints is N(mu, Q^-1)", {
  cm <- cblocks(n, 2, G_list)
  S <- Reduce(`+`, Map(function(t, g) exp(-t) * g, theta, G_list))
  set.seed(1); x <- rnorm(n)
  expect_equal(gmrf_log_density(cm, x, theta),
               -0.5 * n * log(2 * pi) - 0.5 * as.numeric(determinant(S)$modulus) -
                 0.5 * sum(x * solve(S, x)))
  X <- cbind(x, -x)
  expect_equal(gmrf_log_density(cm, X, theta),
               c(gmrf_log_density(cm, x, theta), gmrf_log_density(cm, -x, theta)))
})

test_that("log density and samples with constraints, Q_add and b", {
  cm <- cblocks(n, 2, G_list)
  cx <- list(A = rbind(rep(1, n), c(1, -1, rep(0, n - 2))), e = c(0.5, 0))
  Q <- gmrf_Q(cm, theta)
  set.seed(2)
  x <- gmrf_sample(cm, theta, n = 3, constr = cx)
  expect_equal(as.numeric(cx$A %*% x), rep(cx$e, 3), tolerance = 1e-10)
  expect_equal(gmrf_log_density(cm, x[, 1], theta, constr = cx),
               ld_sub(Q, numeric(n), cx$A, cx$e, x[, 1]), tolerance = 1e-8)
  expect_equal(gmrf_log_density(cm, x[, 1] + 1e-3, theta, constr = cx), -Inf)

  Q_add <- diag(4, n); b <- 4 * seq_len(n) / n
  Qp <- as.matrix(Q) + Q_add
  expect_equal(gmrf_log_density(cm, x[, 1], theta, constr = cx, Q_add = Q_add, b = b),
               ld_sub(Qp, solve(Qp, b), cx$A, cx$e, x[, 1]), tolerance = 1e-8)
  ## z = 0 gives the constrained mean
  s <- gmrf_sample(cm, theta, constr = cx, z = numeric(n), Q_add = Q_add, b = b)[, 1]
  V <- solve(Qp); m <- solve(Qp, b)
  ref <- m - V %*% t(cx$A) %*% solve(cx$A %*% V %*% t(cx$A), cx$A %*% m - cx$e)
  expect_equal(s, as.numeric(ref), tolerance = 1e-8)
})

test_that("log prior: the density by default, INLA's value (+ K) with inla = TRUE", {
  skip_if_not_installed("INLA")
  W <- matrix(0, 6, 6)
  W[cbind(1:5, 2:6)] <- 1
  W <- W + t(W)
  m <- fbesag(W, id = rep(1, 6), nrep = 2)
  lam <- -log(1e-5)
  th <- 0.7
  expect_equal(gmrf_log_prior(m, th), log(lam / 2) - 0.5 * th - lam * exp(-0.5 * th))
  expect_equal(gmrf_log_prior(m, th, inla = TRUE) - gmrf_log_prior(m, th),
               m$f$cgeneric$data$doubles$K)
})

test_that("intrinsic models: samples on the constraints, density on the subspace", {
  skip_if_not_installed("INLA")
  op <- options(mc.cores = 1)
  on.exit(options(op))
  W <- matrix(0, 6, 6)
  W[cbind(1:5, 2:6)] <- 1
  W <- W + t(W)
  for (m in list(fbesag(W, id = c(1, 1, 1, 2, 2, 2)), fbym2(W, id = c(1, 1, 1, 2, 2, 2)))) {
    th <- gmrf_initial(m) + 0.2
    cx <- m$f$extraconstr
    set.seed(3)
    x <- gmrf_sample(m, th, n = 2)
    expect_equal(as.numeric(cx$A %*% x), rep(cx$e, 2), tolerance = 1e-9)
    expect_equal(gmrf_log_density(m, x[, 1], th),
                 ld_sub(gmrf_Q(m, th), gmrf_mu(m, th), cx$A, cx$e, x[, 1]), tolerance = 1e-7)
  }
})
