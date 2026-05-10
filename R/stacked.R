#' Stacked block model
#'
#' Latent dimension is `p * n`: blocks are concatenated.  The precision
#' is block-diagonal `diag(exp(theta_1) R_1, ..., exp(theta_p) R_p)`.
#' Each `R_j` is supplied as a dense `n x n` precision/structure matrix
#' (typically `R_j = G_j^{-1}` for some block covariance `G_j`).
#'
#' @param n block size.
#' @param p number of blocks.
#' @param R_list a length-`p` list of `n x n` per-block precision matrices.
#' @param tau0 horseshoe-style global scale (default 1). See
#'   [cblocks()].
#' @param shlib optional path to a specific `.so`/`.dll`. When supplied
#'   it overrides the resolution chain (local dir, INLA, package).
#' @param prefer_inla if TRUE, prefer INLA's precompiled library when
#'   shipped.
#' @param debug logical.
#' @return cgeneric model handle.
#' @export
sblocks <- function(n, p, R_list, tau0 = 1.0, shlib = NULL,
                          prefer_inla = TRUE, debug = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("INLA is required", call. = FALSE)
  n <- as.integer(n); p <- as.integer(p)
  if (!is.list(R_list) || length(R_list) != p)
    stop(sprintf("R_list must be a list of length p = %d", p), call. = FALSE)
  if (!is.numeric(tau0) || length(tau0) != 1L || tau0 <= 0)
    stop("tau0 must be a positive scalar", call. = FALSE)

  ## Stack the per-block precisions into one (p*n) x n dense matrix,
  ## row-major (each block j's R_j occupies rows (j-1)*n+1..j*n).
  ## This matches the layout collapsed.c / partial.c expect for G.
  R_stacked  <- matrix(0.0, p * n, n)
  logdet_R_all <- 0.0
  for (j in seq_len(p)) {
    Rj <- as.matrix(R_list[[j]])
    if (!identical(dim(Rj), c(n, n)))
      stop(sprintf("R_list[[%d]] must be %d x %d", j, n, n), call. = FALSE)
    ch <- tryCatch(chol(Rj), error = function(e) NULL)
    if (is.null(ch))
      stop(sprintf("R_list[[%d]] must be positive-definite", j), call. = FALSE)
    R_stacked[((j - 1L) * n + 1L):(j * n), ] <- Rj
    logdet_R_all <- logdet_R_all + 2 * sum(log(diag(ch)))
  }

  INLA::inla.cgeneric.define(
    model        = "inla_cgeneric_sblocks",
    shlib        = gmrfs_resolve_shlib("sblocks", shlib, prefer_inla),
    n            = as.integer(n * p),
    p            = p,
    R_stacked    = R_stacked,
    logdet_R_all = as.numeric(logdet_R_all),
    tau0         = as.numeric(tau0),
    debug        = isTRUE(debug)
  )
}
