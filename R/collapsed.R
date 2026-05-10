#' Collapsed block model
#'
#' Latent dimension stays at `n` regardless of `p`: INLA sees a single
#' `f(idx, model = cm)` instead of a stacked p-fold chain. The covariance
#' is `Sigma(theta) = sum_j exp(-theta_j) G_j`.
#'
#' @param n block size (= number of individuals; latent dim).
#' @param p number of blocks (= length of `G_list`).
#' @param G_list a length-`p` list of `n x n` per-block covariance matrices.
#' @param tau0 horseshoe-style global scale on the standard-deviation scale of
#'   each block (default 1). Larger tau0 -> looser prior on sigma_j; smaller
#'   tau0 -> more aggressive shrinkage. The C-side prior reads tau0 from
#'   `data->doubles` and falls back to 1 when not set.
#' @param shlib optional path to a specific `.so`/`.dll`. When supplied it
#'   overrides the resolution chain (local dir, INLA, package).
#' @param prefer_inla if TRUE, prefer INLA's precompiled library when shipped.
#' @param debug logical.
#' @return cgeneric model handle.  After fitting, `theta[j]` corresponds to
#'   `G_list[[j]]` (no reordering).
#' @export
cblocks <- function(n, p, G_list, tau0 = 1.0, shlib = NULL,
                    prefer_inla = TRUE, debug = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("INLA is required", call. = FALSE)
  n <- as.integer(n); p <- as.integer(p)
  if (!is.list(G_list) || length(G_list) != p)
    stop(sprintf("G_list must be a list of length p = %d", p), call. = FALSE)
  if (!is.numeric(tau0) || length(tau0) != 1L || tau0 <= 0)
    stop("tau0 must be a positive scalar", call. = FALSE)

  ## Stack into one (p*n) x n dense matrix; row block j is G_list[[j]].
  G <- matrix(0.0, p * n, n)
  for (j in seq_len(p)) {
    Gj <- as.matrix(G_list[[j]])
    if (!identical(dim(Gj), c(n, n)))
      stop(sprintf("G_list[[%d]] must be %d x %d", j, n, n), call. = FALSE)
    G[((j - 1L) * n + 1L):(j * n), ] <- Gj
  }

  INLA::inla.cgeneric.define(
    model = "inla_cgeneric_collapsed",
    shlib = gmrfs_resolve_shlib("collapsed", shlib, prefer_inla),
    n     = n,
    p     = p,
    G     = G,
    tau0  = as.numeric(tau0),
    debug = isTRUE(debug)
  )
}
