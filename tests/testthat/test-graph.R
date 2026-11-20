## graph reading and scaling without INLA (R/graph.R, src/qinv.c)

## triangle {1,2,3}, path {4,5,6}, island 7
W <- matrix(0, 7, 7)
for (e in list(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6))) W[e[1], e[2]] <- W[e[2], e[1]] <- 1

test_that("graphs from a matrix, a file and an nb object agree", {
  g <- gmrfs_read_graph(W)
  expect_equal(g$n, 7)
  expect_equal(g$nnbs, c(2, 2, 2, 1, 2, 1, 0))
  expect_equal(lapply(g$cc$nodes, as.integer), list(1:3, 4:6, 7L))
  f <- tempfile()
  on.exit(unlink(f))
  ## a graph file, neighbours in any order
  writeLines(c("7", "1 2 3 2", "2 2 3 1", "3 2 2 1", "4 1 5", "5 2 6 4", "6 1 5", "7 0"), f)
  gf <- gmrfs_read_graph(f)
  expect_equal(lapply(gf$nbs, function(x) sort(as.integer(x))), lapply(g$nbs, function(x) sort(as.integer(x))))
  nb <- structure(lapply(g$nbs, function(x) if (length(x)) as.integer(x) else 0L), class = "nb")
  gn <- gmrfs_read_graph(nb)
  expect_equal(lapply(gn$nbs, as.integer), lapply(g$nbs, as.integer))
  expect_equal(as.matrix(gmrfs_graph2matrix(g)), (W != 0) + diag(7), ignore_attr = TRUE)
})

test_that("marginal variances by the Takahashi recursions = dense inverse", {
  set.seed(1)
  B <- matrix(rnorm(144), 12)
  Q <- Matrix::Matrix(crossprod(B) * (abs(row(diag(12)) - col(diag(12))) < 3) + diag(12), sparse = TRUE)
  V <- solve(as.matrix(Q))
  expect_equal(gmrfs_qinv_diag(Q), diag(V), tolerance = 1e-10)
  A <- rbind(rep(1, 12))
  Vc <- V - V %*% t(A) %*% solve(A %*% V %*% t(A), A %*% V)
  expect_equal(gmrfs_qinv_diag(Q, list(A = A, e = 0)), diag(Vc), tolerance = 1e-10)
})

test_that("scaling: geometric mean of the constrained marginal variances is 1", {
  R <- diag(rowSums(W)) - W
  S <- gmrfs_scale_model(R, list(A = matrix(1, 1, 7), e = 0))
  for (i in list(1:3, 4:6)) {
    Si <- as.matrix(S[i, i]); n <- length(i); J <- matrix(1 / n, n, n)
    ## covariance under sum-to-zero: (S + J/n)^-1 - J/n, with 1 the null vector of S
    v <- diag(solve(Si + J) - J)
    expect_equal(exp(mean(log(v))), 1, tolerance = 1e-6)
  }
  expect_equal(S[7, 7], 1)
})

test_that("fbesag and fbym2 are built without INLA's functions", {
  op <- options(mc.cores = 1)
  on.exit(options(op))
  fm <- fbesag(W, id = c(1, 1, 2, 2, 2, 1, 2))
  bm <- fbym2(W, id = c(1, 1, 2, 2, 2, 1, 2))
  ## components {1,2,3} and {4,5,6}: one sum-to-zero each, the island none
  expect_equal(nrow(fm$f$extraconstr$A), 2)
  expect_equal(nrow(bm$f$extraconstr$A), 2)
  x <- gmrf_sample(fm, c(0.2, -0.1))
  expect_equal(as.numeric(fm$f$extraconstr$A %*% x), c(0, 0), tolerance = 1e-10)
  expect_true(is.finite(gmrf_log_density(bm, gmrf_sample(bm), gmrf_initial(bm))))
})
