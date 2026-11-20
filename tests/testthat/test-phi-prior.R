## the prior of logit(phi) in fbym2: INLA's table inside its range, the exact
## tail of the PC prior beyond it, and a proper density on the whole line

## triangle {1,2,3}, path {4,5,6}, island 7: rank deficiency 2
W <- matrix(0, 7, 7)
for (e in list(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6))) W[e[1], e[2]] <- W[e[2], e[1]] <- 1

test_that("the prior of logit(phi) is continuous, decreasing beyond the table and proper", {
  op <- options(mc.cores = 1)
  on.exit(options(op))
  m <- fbym2(W, rep(1, 7))
  dd <- m$f$cgeneric$data$doubles
  tmax <- dd$phi_tail[1]
  expect_equal(dd$phi_tail[3], 2)
  lp <- function(t2, inla = TRUE) vapply(t2, function(t) gmrf_log_prior(m, c(0, t), inla = inla), 0)

  ## continuous at the end of the table, and equal to the table inside it
  expect_equal(lp(tmax + 1e-9), lp(tmax - 1e-9), tolerance = 1e-6)
  i <- c(10, 5000, 9990)
  y_tab <- vapply(dd$phi_x[i], function(t) lp(t) - lp(dd$phi_x[i[1]]), 0)
  expect_equal(y_tab, dd$phi_y[i] - dd$phi_y[i[1]], tolerance = 1e-6)

  ## the tail decreases, and has the exact form: d^2 = d_max^2 + r (L - L_max);
  ## the rest of the log prior is pc.prec(theta = 0) and K
  expect_true(all(diff(lp(c(12, 20, 50, 200, 1000))) < 0))
  lpc <- log(-log(0.01) / 2) + log(0.01)
  t <- 30
  L <- function(x) x + log1p(exp(-x))
  d <- sqrt(dd$phi_tail[2]^2 + 2 * (L(t) - L(tmax)))
  expect_equal(lp(t), dd$phi_tail[5] - dd$phi_tail[4] * d + log(2 * plogis(t) / d) + lpc + dd$K,
               tolerance = 1e-8)

  ## the exact prior (inla = FALSE) integrates to 1 over logit(phi):
  ## numerically up to T, analytically beyond (int exp(c0 - lambda d) 2 dd)
  T <- 1e4
  x <- c(seq(-60, tmax, length.out = 20001), exp(seq(log(tmax), log(T), length.out = 20001))[-1])
  y <- exp(lp(x, inla = FALSE) - lpc)
  inside <- sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
  dT <- sqrt(dd$phi_tail[2]^2 + 2 * (L(T) - L(tmax)))
  beyond <- 2 * exp(dd$phi_tail[5] - dd$phi_tail[4] * dT) / dd$phi_tail[4] / exp(dd$phi_log_mass)
  expect_equal(inside + beyond, 1, tolerance = 1e-4)
})

test_that("Q stays finite when phi is close to 1", {
  op <- options(mc.cores = 1)
  on.exit(options(op))
  m <- fbym2(W, rep(1, 7))
  Q <- gmrf_Q(m, c(0, 40))
  expect_true(all(is.finite(Q@x)))
})
