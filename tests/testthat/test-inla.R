## Agreement with INLA; skipped where INLA is not installed (e.g. on CRAN).

ar1_cov <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))

test_that("the local model handle matches inla.cgeneric.define()", {
  skip_if_not_installed("INLA")
  S <- Matrix::sparseMatrix(i = c(1, 1, 2), j = c(1, 2, 2), x = c(2, -1, 3))
  args <- list(model = "inla_cgeneric_gmrfs_cblocks", shlib = gmrfs:::gmrfs_own_shlib(),
               n = 2L, debug = FALSE, p = 1L, G = matrix(c(1, 2, 3, 4, 5, 6), 3, 2),
               S = gmrfs:::gmrfs_as_triplet(S), tau0 = 1.5, v = c(1, 2))
  expect_equal(do.call(gmrfs:::.gmrfs_define_local, args),
               do.call(INLA::inla.cgeneric.define, args))
})

test_that("gmrf_* agree with INLA's own evaluation", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  op <- options(mc.cores = 1)
  on.exit(options(op))
  W <- matrix(0, 6, 6)
  W[cbind(1:5, 2:6)] <- 1
  W <- W + t(W)
  models <- list(
    cblocks = cblocks(6, 2, list(ar1_cov(6, 0.5), ar1_cov(6, 0.2)), prefer_inla = FALSE),
    fbesag  = fbesag(W, id = c(1, 1, 1, 2, 2, 2), prefer_inla = FALSE),
    fbym2   = fbym2(W, id = c(1, 1, 1, 2, 2, 2), prefer_inla = FALSE)
  )
  for (nm in names(models)) {
    m <- models[[nm]]
    ref <- INLA::inla.cgeneric.q(m)
    expect_equal(gmrf_initial(m), ref$theta, tolerance = 1e-6, info = nm)
    expect_equal(as.matrix(gmrf_Q(m)), as.matrix(ref$Q), tolerance = 1e-6,
                 ignore_attr = TRUE, info = nm)
    expect_equal(gmrf_log_prior(m, inla = TRUE), ref$log.prior, tolerance = 1e-6, info = nm)
  }
})
