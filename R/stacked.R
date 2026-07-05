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
sblocks <- function(n, p, R_list,
                    tau         = NULL,
                    tau0        = 1.0,
                    noise_scale = FALSE,
                    prior       = "horseshoe",
                    param       = NULL,
                    fixed       = NULL,
                    shlib       = NULL,
                    prefer_inla = TRUE,
                    debug       = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("INLA is required", call. = FALSE)
  n <- as.integer(n); p <- as.integer(p)
  if (!is.list(R_list) || length(R_list) != p)
    stop(sprintf("R_list must be a list of length p = %d", p), call. = FALSE)

  ## --- resolve the global-scale option `tau` (see cblocks() for logic) ---
  if (is.null(tau)) tau <- tau0
  marg_tau <- 0L
  tau_num  <- 1.0
  if (is.character(tau) && length(tau) == 1L && identical(tau, "marginal")) {
    marg_tau <- 1L
    tau_num  <- 1.0
  } else if (is.numeric(tau) && length(tau) == 1L && is.finite(tau) && tau > 0) {
    marg_tau <- 0L
    tau_num  <- as.numeric(tau)
  } else {
    stop("`tau` must be a positive number or the string \"marginal\"",
         call. = FALSE)
  }

  pp <- gmrfs_prior_param(prior, param, tau_num)
  fixed_flag  <- if (is.null(fixed)) 0L else 1L
  theta_fixed <- if (is.null(fixed)) 0.0 else log(as.numeric(fixed))

  ## If any new feature is on, force gmrfs's own .so.
  needs_new <- (marg_tau == 1L) || isTRUE(noise_scale)
  if (needs_new) prefer_inla <- FALSE

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
    model        = "inla_cgeneric_gmrfs_sblocks",
    shlib        = gmrfs_resolve_shlib("sblocks", shlib, prefer_inla),
    n            = as.integer(n * p),
    p            = p,
    R_stacked    = R_stacked,
    logdet_R_all = as.numeric(logdet_R_all),
    tau0         = as.numeric(tau_num),
    noise_scale  = as.integer(isTRUE(noise_scale)),
    marg_tau     = as.integer(marg_tau),
    prior_id     = gmrfs_prior_code(prior),
    param        = pp,
    fixed        = fixed_flag,
    theta_fixed  = theta_fixed,
    debug        = isTRUE(debug)
  )
}
