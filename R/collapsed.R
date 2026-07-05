#' Collapsed block model
#'
#' Latent dimension stays at `n` regardless of `p`: INLA sees a single
#' `f(idx, model = cm)` instead of a stacked p-fold chain. The covariance
#' is `Sigma(theta) = sum_j exp(-theta_j) G_j`.
#'
#' @param n block size (= number of individuals; latent dim).
#' @param p number of blocks (= length of `G_list`).
#' @param G_list a length-`p` list of `n x n` per-block covariance matrices.
#' @param tau global scale of the horseshoe. Either a positive number to fix
#'   `tau` at that constant, or the string `"marginal"` to integrate `tau`
#'   out of the prior via a 60-point log-tau quadrature (recommended when
#'   you want a data-driven global scale). Defaults to `tau0` if `tau0` is
#'   supplied and `tau` is not, otherwise `1.0`.
#' @param tau0 back-compat alias for `tau` when `tau` is not supplied.
#'   Older gmrfs code that passed `tau0 = X` continues to work.
#' @param noise_scale logical. If TRUE, read the Gaussian family's noise
#'   precision from INLA at each call and multiply the resulting
#'   `sigma_e` into the horseshoe global scale, so the prior becomes
#'   `sigma_j ~ C+(0, tau * sigma_e)` (matches the sigma^2-scaled
#'   horseshoe). Requires an INLA that populates `data->theta_all_*`.
#'   Ignored if the Gaussian family precision is fixed.  Default FALSE.
#' @param nugget diagonal floor added to `Sigma(theta)` before the
#'   Cholesky, for numerical stability when all blocks are strongly
#'   shrunk. Default `0` for backward compatibility; auto-bumped to
#'   `1e-4` if `tau = "marginal"` (marg-tau visits arbitrarily small
#'   tau on the quadrature grid, where a nugget is essential).
#' @param prior prior name (default `"horseshoe"`). Only used when
#'   `tau` is a number and `noise_scale` is FALSE (the "plain" path).
#'   `noise_scale = TRUE` or `tau = "marginal"` forces the horseshoe.
#' @param param optional prior parameters; if NULL, defaults to the
#'   scale set by `tau`.
#' @param fixed optional vector of length 1 or `p` of fixed block
#'   precisions (on the original, non-log scale). When supplied, the
#'   block log-precisions are held at those values and no
#'   hyperparameters are estimated.
#' @param shlib optional path to a specific `.so`/`.dll`.
#' @param prefer_inla if TRUE, prefer INLA's precompiled library when shipped.
#' @param debug logical.
#' @return cgeneric model handle.
#' @export
cblocks <- function(n, p, G_list,
                    tau         = NULL,
                    tau0        = 1.0,
                    noise_scale = FALSE,
                    nugget      = 0.0,
                    prior       = "horseshoe",
                    param       = NULL,
                    fixed       = NULL,
                    shlib       = NULL,
                    prefer_inla = TRUE,
                    debug       = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("INLA is required", call. = FALSE)
  n <- as.integer(n); p <- as.integer(p)
  if (!is.list(G_list) || length(G_list) != p)
    stop(sprintf("G_list must be a list of length p = %d", p), call. = FALSE)

  ## --- resolve the global-scale option `tau` ---
  ## Priority: explicit tau > tau0 (back-compat) > default 1.0
  if (is.null(tau)) tau <- tau0
  marg_tau <- 0L
  tau_num  <- 1.0
  if (is.character(tau) && length(tau) == 1L && identical(tau, "marginal")) {
    marg_tau <- 1L
    tau_num  <- 1.0                          ## unused in marg path; keep sensible
    if (nugget <= 0.0) nugget <- 1e-4       ## marg-tau needs a stability floor
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

  ## If the user opts into any feature that INLA's precompiled libcblocks
  ## does not have (marg_tau, noise_scale, nugget), fall back to gmrfs's own
  ## shared library so the new code paths are actually executed.
  needs_new <- (marg_tau == 1L) || isTRUE(noise_scale) || (nugget > 0.0)
  if (needs_new) prefer_inla <- FALSE

  ## Stack into one (p*n) x n dense matrix; row block j is G_list[[j]].
  G <- matrix(0.0, p * n, n)
  for (j in seq_len(p)) {
    Gj <- as.matrix(G_list[[j]])
    if (!identical(dim(Gj), c(n, n)))
      stop(sprintf("G_list[[%d]] must be %d x %d", j, n, n), call. = FALSE)
    G[((j - 1L) * n + 1L):(j * n), ] <- Gj
  }

  INLA::inla.cgeneric.define(
    model       = "inla_cgeneric_gmrfs_cblocks",
    shlib       = gmrfs_resolve_shlib("cblocks", shlib, prefer_inla),
    n           = n,
    p           = p,
    G           = G,
    tau0        = as.numeric(tau_num),
    noise_scale = as.integer(isTRUE(noise_scale)),
    marg_tau    = as.integer(marg_tau),
    nugget      = as.numeric(nugget),
    prior_id    = gmrfs_prior_code(prior),
    param       = pp,
    fixed       = fixed_flag,
    theta_fixed = theta_fixed,
    debug       = isTRUE(debug)
  )
}
