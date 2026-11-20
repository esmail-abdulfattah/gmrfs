## fbesag and fbym2 carry diagonal = 1e-8 for INLA's f(): INLA's default of
## 1e-4 enters the normalising constant it computes for them

W <- matrix(0, 5, 5)
W[cbind(1:4, 2:5)] <- 1
W <- W + t(W)

test_that("fbesag and fbym2 set a small diagonal in their handle", {
  op <- options(mc.cores = 1)
  on.exit(options(op))
  expect_equal(fbym2(W, rep(1, 5))$f$diagonal, 1e-8)
  expect_equal(fbesag(W, c(1, 1, 1, 2, 2))$f$diagonal, 1e-8)
})

test_that("the jitter of the intrinsic models is 1e-8", {
  op <- options(mc.cores = 1)
  on.exit(options(op))
  ev <- eigen(as.matrix(gmrf_Q(fbym2(W, rep(1, 5)), c(0, 0))), symmetric = TRUE)$values
  expect_equal(min(ev), 1e-8, tolerance = 1e-6)
})
