## Model handles and their evaluation from R, without INLA.

## Build the model handle. With INLA installed this is
## INLA::inla.cgeneric.define(); without it, the same object is built here
## (same layout and argument order, which the C code relies on), so the
## handle can still be evaluated with gmrf_Q() and friends.
gmrfs_cgeneric_define <- function(model, name, n, ..., shlib = NULL,
                                  prefer_inla = TRUE, debug = FALSE) {
  n <- as.integer(n)
  if (requireNamespace("INLA", quietly = TRUE)) {
    return(INLA::inla.cgeneric.define(
      model = model, shlib = gmrfs_resolve_shlib(name, shlib, prefer_inla, model),
      n = n, debug = debug, ...))
  }

  .gmrfs_define_local(model, gmrfs_resolve_shlib(name, shlib, prefer_inla = FALSE),
                      n, debug, ...)
}

.gmrfs_define_local <- function(model, shlib, n, debug, ...) {
  a <- c(n = n, model = model, shlib = shlib, debug = as.integer(debug), list(...))
  if (any(names(a) == "")) stop("unnamed arguments are not allowed", call. = FALSE)
  data <- list(ints = list(), doubles = list(), characters = list(),
               matrices = list(), smatrices = list())
  for (idx in seq_along(a)) {
    xx <- a[[idx]]
    if (is.matrix(xx)) {
      ## dense matrices are stored row by row
      data$matrices[[names(a)[idx]]] <- c(nrow(xx), ncol(xx), as.numeric(t(xx)))
    } else if (methods::is(xx, "Matrix")) {
      xx <- gmrfs_as_triplet(xx)
      data$smatrices[[names(a)[idx]]] <- c(nrow(xx), ncol(xx), length(xx@x), xx@i, xx@j, xx@x)
    } else if (is.integer(xx)) {
      data$ints <- c(data$ints, a[idx])
    } else if (is.character(xx)) {
      data$characters <- c(data$characters, a[idx])
    } else if (is.double(xx)) {
      data$doubles <- c(data$doubles, a[idx])
    } else {
      stop(sprintf("unknown type of argument '%s'", names(a)[idx]), call. = FALSE)
    }
  }
  structure(
    list(f = list(model = "cgeneric", n = n,
                  cgeneric = structure(list(model = model, shlib = shlib, n = n,
                                            debug = debug, data = data),
                                       class = "inla.cgeneric.f"))),
    class = "inla.cgeneric")
}

## general triplet form (dgTMatrix), all entries stored explicitly
gmrfs_as_triplet <- function(x) {
  methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "TsparseMatrix")
}

## cgeneric commands, as in inla_cgeneric_cmd_tp
.gmrfs_cmd <- c(Q = 1L, graph = 2L, mu = 3L, initial = 4L,
                log_norm_const = 5L, log_prior = 6L)

.gmrfs_eval <- function(model, cmd, theta = numeric(0)) {
  cg <- model$f$cgeneric
  if (!inherits(cg, "inla.cgeneric.f"))
    stop("`model` must be a model handle created by this package", call. = FALSE)
  .Call(gmrfs_cgeneric_eval, cg$model, .gmrfs_cmd[[cmd]], as.numeric(theta), cg$data)
}

## theta defaults to the initial values; its length is checked here because
## the C code reads it without bounds checks
.gmrfs_theta <- function(model, theta) {
  init <- gmrf_initial(model)
  if (is.null(theta)) return(init)
  theta <- as.numeric(theta)
  if (length(theta) != length(init))
    stop(sprintf("`theta` must have length %d", length(init)), call. = FALSE)
  if (anyNA(theta)) stop("`theta` must not contain missing values", call. = FALSE)
  theta
}

#' Evaluate a latent Gaussian Markov random field
#'
#' The models in this package are latent Gaussian fields
#' `x | theta ~ N(mu(theta), Q(theta)^-1)` with a prior on the
#' hyperparameters `theta`. These functions evaluate each ingredient of a
#' model handle directly, so the models can be used in your own code (an
#' MCMC sampler, an optimiser, a simulation) and not only inside INLA.
#' INLA does not need to be installed.
#'
#' `theta` is on the internal scale of the model (log precisions and so
#' on; see the help page of each model). Linear constraints that belong to
#' a model, if any, are in `model$f$extraconstr`.
#'
#' The evaluation always uses the code compiled with this package, also
#' when the handle points to a library shipped by INLA. Options that need
#' information from a running INLA fit (`noise_scale = TRUE`) are ignored.
#'
#' @param model a model handle, as returned by [cblocks()], [pblocks()],
#'   [sblocks()], [fbesag()] or [fbym2()].
#' @param theta hyperparameters on the internal scale. Defaults to the
#'   initial values, `gmrf_initial(model)`.
#' @param inla for `gmrf_log_prior()`: if `FALSE` (default), the log prior
#'   density of `theta`. If `TRUE`, the value the C code returns to INLA:
#'   for [fbesag()] and [fbym2()] this adds a constant `K` that puts INLA's
#'   marginal likelihood on the convention of `besag`/`bym2` (INLA leaves
#'   some terms of their normalising constant out); `K` does not depend on
#'   `theta`. For [fbym2()] it also leaves the prior of `logit(phi)`
#'   normalised on its table, as INLA's `bym2`, rather than on the whole
#'   line. The block models have no such constants.
#' @return
#' * `gmrf_Q()`: the precision matrix, a sparse symmetric matrix
#'   (class `dsCMatrix`).
#' * `gmrf_graph()`: the sparsity pattern of the precision matrix, a sparse
#'   symmetric pattern matrix.
#' * `gmrf_mu()`: the mean, a numeric vector of the latent dimension.
#' * `gmrf_initial()`: the initial values of `theta`; its length is the
#'   number of hyperparameters.
#' * `gmrf_log_prior()`: the log prior density of `theta` (with
#'   `inla = TRUE`, INLA's value; see `inla`).
#' * `gmrf_log_norm_const()`: the log normalising constant of
#'   `x | theta`, or `NA` if the model leaves it to be computed from
#'   `Q(theta)`.
#' @examples
#' ar1_cov <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
#' G_list <- list(ar1_cov(10, 0.5), ar1_cov(10, 0.2))
#' cm <- cblocks(n = 10, p = 2, G_list = G_list)
#'
#' theta <- c(0, 1)
#' Q <- gmrf_Q(cm, theta)
#' ## the covariance is sum_j exp(-theta_j) G_j
#' range(solve(Q) - (exp(-theta[1]) * G_list[[1]] + exp(-theta[2]) * G_list[[2]]))
#'
#' gmrf_initial(cm)
#' gmrf_log_prior(cm, theta)
#' gmrf_log_norm_const(cm, theta)
#'
#' ## a draw from the field
#' x <- solve(chol(Q), rnorm(10))
#' @name gmrf_eval
NULL

#' @rdname gmrf_eval
#' @export
gmrf_Q <- function(model, theta = NULL) {
  ret <- .gmrfs_eval(model, "Q", .gmrfs_theta(model, theta))
  M <- ret[2]
  if (ret[1] == -1) {
    ## values only, in the order of the graph
    g <- .gmrfs_eval(model, "graph")
    n <- g[1]
    i <- g[2 + seq_len(M)]; j <- g[2 + M + seq_len(M)]; x <- ret[2 + seq_len(M)]
  } else {
    n <- ret[1]
    i <- ret[2 + seq_len(M)]; j <- ret[2 + M + seq_len(M)]; x <- ret[2 + 2 * M + seq_len(M)]
  }
  Matrix::sparseMatrix(i = pmin(i, j) + 1, j = pmax(i, j) + 1, x = x,
                       dims = c(n, n), symmetric = TRUE)
}

#' @rdname gmrf_eval
#' @export
gmrf_graph <- function(model) {
  g <- .gmrfs_eval(model, "graph")
  n <- g[1]; M <- g[2]
  i <- g[2 + seq_len(M)]; j <- g[2 + M + seq_len(M)]
  Matrix::sparseMatrix(i = pmin(i, j) + 1, j = pmax(i, j) + 1,
                       dims = c(n, n), symmetric = TRUE)
}

#' @rdname gmrf_eval
#' @export
gmrf_mu <- function(model, theta = NULL) {
  ret <- .gmrfs_eval(model, "mu", .gmrfs_theta(model, theta))
  ## a leading 0 means a zero mean
  if (ret[1] == 0) numeric(model$f$n) else ret[-1]
}

#' @rdname gmrf_eval
#' @export
gmrf_initial <- function(model) {
  .gmrfs_eval(model, "initial")[-1]
}

#' @rdname gmrf_eval
#' @export
gmrf_log_prior <- function(model, theta = NULL, inla = FALSE) {
  ret <- .gmrfs_eval(model, "log_prior", .gmrfs_theta(model, theta))
  if (is.null(ret)) return(NA_real_)
  ## fbesag/fbym2 pass K (already times nrep) to the C code, which adds it
  K <- model$f$cgeneric$data$doubles$K
  if (!isTRUE(inla) && !is.null(K)) ret <- ret - K
  ## fbym2: the tabulated prior of logit(phi) integrates to 1 on the table,
  ## as INLA's; divide by its mass on the whole line for the exact prior
  log_mass <- model$f$cgeneric$data$doubles$phi_log_mass
  if (!isTRUE(inla) && !is.null(log_mass)) ret <- ret - log_mass
  ret
}

#' @rdname gmrf_eval
#' @export
gmrf_log_norm_const <- function(model, theta = NULL) {
  ret <- .gmrfs_eval(model, "log_norm_const", .gmrfs_theta(model, theta))
  if (is.null(ret)) NA_real_ else ret
}

## Sparse Cholesky of Q(theta) (+ Q_add) and, with constraints A x = e, the
## pieces of the distribution of A x: Q^-1 A' and S = A Q^-1 A'. Shared by
## gmrf_log_density() and gmrf_sample().
.gmrfs_factor <- function(model, theta, constr, Q_add = NULL, b = NULL) {
  Q  <- gmrf_Q(model, theta)
  mu <- gmrf_mu(model, theta)
  if (!is.null(Q_add) || !is.null(b)) {
    ## x | data ~ N((Q + Q_add)^-1 (Q mu + b), (Q + Q_add)^-1)
    b <- if (is.null(b)) numeric(nrow(Q)) else as.numeric(b)
    if (length(b) != nrow(Q)) stop(sprintf("`b` must have length %d", nrow(Q)), call. = FALSE)
    Qmu_b <- as.numeric(Q %*% mu) + b
    if (!is.null(Q_add)) {
      if (!all(dim(Q_add) == dim(Q)))
        stop(sprintf("`Q_add` must be %d x %d", nrow(Q), nrow(Q)), call. = FALSE)
      Q <- Matrix::forceSymmetric(Q + Matrix::Matrix(Q_add, sparse = TRUE), uplo = "U")
    }
    ch <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE, super = FALSE)
    mu <- as.numeric(Matrix::solve(ch, Qmu_b))
  } else {
    ch <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE, super = FALSE)
  }
  out <- list(Q = Q, mu = mu, ch = ch, A = NULL)
  if (!is.null(constr) && NROW(constr$A) > 0) {
    A <- matrix(as.numeric(constr$A), ncol = nrow(Q))
    e <- as.numeric(constr$e)
    if (length(e) != nrow(A))
      stop("`constr$e` must have one value per row of `constr$A`", call. = FALSE)
    QiAt <- as.matrix(Matrix::solve(out$ch, t(A)))
    out$A <- A; out$e <- e; out$QiAt <- QiAt
    out$S <- A %*% QiAt
  }
  out
}

## x as an n x m matrix, one column per vector
.gmrfs_as_cols <- function(x, n, what = "x") {
  x <- as.matrix(x)
  if (nrow(x) != n) {
    if (ncol(x) == n && nrow(x) == 1L) x <- t(x)
    else stop(sprintf("`%s` must have length %d (or %d rows)", what, n, n), call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

#' Log density and samples of the latent field
#'
#' `gmrf_log_density()` evaluates `log p(x | theta)` and `gmrf_sample()`
#' draws `x | theta`, for any model in this package, including its linear
#' constraints `A x = e`. With [gmrf_log_prior()] this is all an MCMC
#' sampler needs (see `demo("fbym2_mcmc", package = "gmrfs")`). INLA does not
#' need to be installed.
#'
#' Both use the precision matrix `Q` of [gmrf_Q()] and a sparse Cholesky
#' factorisation. Without constraints, `x | theta ~ N(mu, Q^-1)`. With
#' constraints, `x` is Gaussian on the subspace `A x = e`, as in INLA: the
#' density, with respect to Lebesgue measure on that subspace, is
#' `pi(x) / pi(A x) * |A A'|^(-1/2)`, where `pi(x)` is the `N(mu, Q^-1)`
#' density and `pi(A x)` that of `A x`; samples are drawn from
#' `N(mu, Q^-1)` and then corrected to `A x = e` (conditioning by kriging).
#'
#' With `Q_add` and/or `b`, both functions work instead with
#' `N((Q + Q_add)^-1 (Q mu + b), (Q + Q_add)^-1)` under the same constraints:
#' the distribution of `x | theta, y` for Gaussian data `y`. With
#' `y ~ N(B x, Q_y^-1)`, `Q_add = t(B) %*% Q_y %*% B` and
#' `b = t(B) %*% Q_y %*% y`. This gives exact Gibbs steps for `x`, and, at
#' any `x` on the constraints, the marginal likelihood
#' `log p(y | theta) = gmrf_log_density(model, x, theta) + log p(y | x) -
#' gmrf_log_density(model, x, theta, Q_add = Q_add, b = b)`.
#'
#' The intrinsic models ([fbesag()], [fbym2()]) add a small jitter (`1e-8`)
#' to the diagonal of `Q`, as INLA does; with the default constraints, which
#' remove the null space of `Q`, the result is the exact intrinsic density
#' up to a relative error of the order of the jitter. The density is
#' complete: it does not leave out the constants that INLA leaves out of
#' `besag`/`bym2`.
#'
#' @inheritParams gmrf_eval
#' @param x the latent vector, or a matrix with one vector per column.
#' @param constr linear constraints `list(A, e)`; by default those of the
#'   model, `model$f$extraconstr` (none for the block models). `NULL` drops
#'   them.
#' @param tol a vector with `max |A x - e|` above `tol * (1 + max |x|)` is
#'   outside the support, and its log density is `-Inf`.
#' @param Q_add optional matrix added to the precision (for example
#'   `t(B) %*% Q_y %*% B` from Gaussian data); see Details.
#' @param b optional vector added to the canonical mean `Q mu`; see
#'   Details.
#' @param n number of samples.
#' @param z optional standard normal vectors to transform, one per column
#'   (then `n` is ignored). For a fixed `z` the sample is a deterministic,
#'   smooth function of `theta`: a non-centred parametrisation, which often
#'   mixes much better in MCMC than updating `theta` given `x`.
#' @return `gmrf_log_density()`: the log density, one value per vector.
#'   `gmrf_sample()`: a matrix with one sample per column.
#' @examples
#' ar1_cov <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
#' cm <- cblocks(n = 10, p = 2, G_list = list(ar1_cov(10, 0.5), ar1_cov(10, 0.2)))
#' theta <- c(0, 1)
#' x <- gmrf_sample(cm, theta, n = 3)
#' gmrf_log_density(cm, x, theta)
#' @examplesIf requireNamespace("INLA", quietly = TRUE)
#' ## an intrinsic model: the samples satisfy the sum-to-zero constraint
#' W <- matrix(0, 5, 5)
#' W[cbind(1:4, 2:5)] <- 1
#' W <- W + t(W)
#' fm <- fbesag(W, id = c(1, 1, 1, 2, 2))
#' x <- gmrf_sample(fm, c(0, 1), n = 2)
#' colSums(x)
#' gmrf_log_density(fm, x, c(0, 1))
#' @name gmrf_density
NULL

#' @rdname gmrf_density
#' @export
gmrf_log_density <- function(model, x, theta = NULL, constr = model$f$extraconstr,
                             tol = 1e-8, Q_add = NULL, b = NULL) {
  f <- .gmrfs_factor(model, theta, constr, Q_add, b)
  N <- length(f$mu)
  x <- unname(.gmrfs_as_cols(x, N))
  r <- x - f$mu
  ld <- -0.5 * N * log(2 * pi) +
    0.5 * as.numeric(Matrix::determinant(f$ch, logarithm = TRUE, sqrt = FALSE)$modulus) -
    0.5 * colSums(r * as.matrix(f$Q %*% r))
  if (is.null(f$A)) return(ld)

  ## subtract log pi(A x) and 0.5 log |A A'|
  k <- nrow(f$A)
  Ax <- f$A %*% x
  cS <- chol(f$S)
  z <- backsolve(cS, Ax - as.numeric(f$A %*% f$mu), transpose = TRUE)
  ld_A <- -0.5 * k * log(2 * pi) - sum(log(diag(cS))) - 0.5 * colSums(z * z)
  ld <- ld - ld_A - 0.5 * as.numeric(determinant(tcrossprod(f$A), logarithm = TRUE)$modulus)
  off <- apply(abs(Ax - f$e), 2, max) > tol * (1 + apply(abs(x), 2, max))
  ld[off] <- -Inf
  ld
}

#' @rdname gmrf_density
#' @export
gmrf_sample <- function(model, theta = NULL, n = 1L, constr = model$f$extraconstr,
                        z = NULL, Q_add = NULL, b = NULL) {
  f <- .gmrfs_factor(model, theta, constr, Q_add, b)
  N <- length(f$mu)
  ## P Q P' = L L', so x = P' L^-T z has covariance Q^-1
  z <- if (is.null(z)) matrix(stats::rnorm(N * n), N, n) else .gmrfs_as_cols(z, N, "z")
  x <- as.matrix(Matrix::solve(f$ch, Matrix::solve(f$ch, z, system = "Lt"), system = "Pt")) + f$mu
  if (!is.null(f$A)) {
    x <- x - f$QiAt %*% solve(f$S, f$A %*% x - f$e)
    ## a second pass removes the rounding error of the first: with the small
    ## jitter of the intrinsic models the draw is large along the null space
    x <- x - f$QiAt %*% solve(f$S, f$A %*% x - f$e)
  }
  dimnames(x) <- NULL
  x
}
