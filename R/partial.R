#' Partial block model: one block kept separate, the rest collapsed
#'
#' Latent dimension is `2 * n`. Block `separate` is treated individually
#' with precision `exp(theta[1]) * R_separate`; the remaining `p - 1`
#' blocks are collapsed into `Sigma_rest = sum_{j != separate}
#' exp(-theta[j]) G_j`. There are `p` hyperparameters: `theta[1]` is the
#' separated block, `theta[2..p]` are the collapsed blocks in their
#' original order minus the separated one.
#'
#' @param n block size.
#' @param p number of blocks (= length of `G_list`).
#' @param G_list a length-`p` list of `n x n` per-block covariance matrices.
#' @param separate integer in `1..p` -- which block to keep separate
#'   (default 1).
#' @inheritParams cblocks
#' @param shlib optional path to a specific `.so`/`.dll`. When supplied it
#'   overrides the resolution chain (local dir, INLA, package).
#' @param prefer_inla if TRUE, prefer INLA's precompiled library when shipped.
#' @param debug logical.
#' @return A model handle, as for [cblocks()]. The attribute `block_order` records
#'   how `theta[j]` maps back to the original `G_list` index --
#'   `theta[j]` corresponds to `G_list[[ block_order[j] ]]`.
#' @examples
#' ar1_cov <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
#' G_list <- list(ar1_cov(10, 0.5), ar1_cov(10, 0.2))
#' pm <- pblocks(n = 10, p = 2, G_list = G_list, separate = 2)
#' attr(pm, "block_order")
#' dim(gmrf_Q(pm, theta = c(0, 1)))
#' @export
pblocks <- function(n, p, G_list,
                    separate    = 1L,
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
  n <- as.integer(n); p <- as.integer(p)
  separate <- as.integer(separate)
  if (length(separate) != 1L || separate < 1L || separate > p)
    stop(sprintf("`separate` must be an integer in 1..%d", p), call. = FALSE)
  if (!is.list(G_list) || length(G_list) != p)
    stop(sprintf("G_list must be a list of length p = %d", p), call. = FALSE)

  ## --- resolve the global-scale option `tau` (see cblocks() for logic) ---
  if (is.null(tau)) tau <- tau0
  marg_tau <- 0L
  tau_num  <- 1.0
  if (is.character(tau) && length(tau) == 1L && identical(tau, "marginal")) {
    marg_tau <- 1L
    tau_num  <- 1.0
    if (nugget <= 0.0) nugget <- 1e-4
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

  ## If any new feature is on, force gmrfs's own .so (INLA's built-in
  ## copy of the old symbol name would ignore these options).
  needs_new <- (marg_tau == 1L) || isTRUE(noise_scale) || (nugget > 0.0)
  if (needs_new) prefer_inla <- FALSE

  ## Reorder so the chosen block goes first; remember the mapping.
  block_order <- c(separate, setdiff(seq_len(p), separate))
  G_list      <- G_list[block_order]

  ## Stack into one (p*n) x n dense matrix in the new order.
  G <- matrix(0.0, p * n, n)
  for (j in seq_len(p)) {
    Gj <- as.matrix(G_list[[j]])
    if (!identical(dim(Gj), c(n, n)))
      stop(sprintf("G_list[[%d]] must be %d x %d (after reordering)",
                   j, n, n), call. = FALSE)
    G[((j - 1L) * n + 1L):(j * n), ] <- Gj
  }

  ## R_separate (= G_list[[separate]]^{-1}) is needed both as the sparse
  ## precision for the top-left block and as its log-determinant.
  R1_dense  <- solve(as.matrix(G_list[[1]]))
  R1_chol   <- tryCatch(chol(R1_dense), error = function(e) NULL)
  if (is.null(R1_chol))
    stop("Separated block must be positive-definite", call. = FALSE)
  R1_logdet <- 2 * sum(log(diag(R1_chol)))

  R1 <- gmrfs_as_triplet(R1_dense)
  ## upper-tri + row-major, as the C side expects.
  keep <- R1@i <= R1@j
  R1@i <- R1@i[keep]; R1@j <- R1@j[keep]; R1@x <- R1@x[keep]
  ord  <- order(R1@i, R1@j)
  R1@i <- R1@i[ord];  R1@j <- R1@j[ord];  R1@x <- R1@x[ord]

  cm <- gmrfs_cgeneric_define(
    "inla_cgeneric_gmrfs_pblocks", "pblocks",
    n           = 2L * n,
    p           = p,
    G           = G,
    R1          = R1,
    R1_logdet   = as.numeric(R1_logdet),
    tau0        = as.numeric(tau_num),
    noise_scale = as.integer(isTRUE(noise_scale)),
    marg_tau    = as.integer(marg_tau),
    nugget      = as.numeric(nugget),
    prior_id    = gmrfs_prior_code(prior),
    param       = pp,
    fixed       = fixed_flag,
    theta_fixed = theta_fixed,
    shlib       = shlib,
    prefer_inla = prefer_inla,
    debug       = isTRUE(debug)
  )
  attr(cm, "block_order") <- block_order
  cm
}
