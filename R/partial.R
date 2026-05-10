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
#' @param tau0 horseshoe-style global scale (default 1). See [cblocks()].
#' @param shlib optional path to a specific `.so`/`.dll`. When supplied it
#'   overrides the resolution chain (local dir, INLA, package).
#' @param prefer_inla if TRUE, prefer INLA's precompiled library when shipped.
#' @param debug logical.
#' @return cgeneric model handle. The attribute `block_order` records
#'   how `theta[j]` maps back to the original `G_list` index --
#'   `theta[j]` corresponds to `G_list[[ block_order[j] ]]`.
#' @export
pblocks <- function(n, p, G_list, separate = 1L, tau0 = 1.0,
                    shlib = NULL, prefer_inla = TRUE, debug = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("INLA is required", call. = FALSE)
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Matrix is required", call. = FALSE)
  n <- as.integer(n); p <- as.integer(p)
  separate <- as.integer(separate)
  if (length(separate) != 1L || separate < 1L || separate > p)
    stop(sprintf("`separate` must be an integer in 1..%d", p), call. = FALSE)
  if (!is.list(G_list) || length(G_list) != p)
    stop(sprintf("G_list must be a list of length p = %d", p), call. = FALSE)
  if (!is.numeric(tau0) || length(tau0) != 1L || tau0 <= 0)
    stop("tau0 must be a positive scalar", call. = FALSE)

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

  R1 <- INLA::inla.as.sparse(R1_dense)
  ## upper-tri + row-major, as the C side expects.
  keep <- R1@i <= R1@j
  R1@i <- R1@i[keep]; R1@j <- R1@j[keep]; R1@x <- R1@x[keep]
  ord  <- order(R1@i, R1@j)
  R1@i <- R1@i[ord];  R1@j <- R1@j[ord];  R1@x <- R1@x[ord]

  cm <- INLA::inla.cgeneric.define(
    model     = "inla_cgeneric_partial",
    shlib     = gmrfs_resolve_shlib("partial", shlib, prefer_inla),
    n         = 2L * n,
    p         = p,
    G         = G,
    R1        = R1,
    R1_logdet = as.numeric(R1_logdet),
    tau0      = as.numeric(tau0),
    debug     = isTRUE(debug)
  )
  attr(cm, "block_order") <- block_order
  cm
}
