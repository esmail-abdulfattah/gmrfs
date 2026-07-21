#' Select candidate regions by Orthogonal Matching Pursuit (OMP)
#'
#' Greedy, INLA-free selection of the regions most likely to carry signal,
#' intended as the cheap first stage before an exact [cblocks()] / [pblocks()]
#' refit on just the selected regions.
#'
#' At each step OMP scores every not-yet-selected region by how much of the
#' current residual lies in its leading eigenspace, selects the best region,
#' \emph{removes} that region's contribution from the residual, and repeats.
#' Removing each selected region before scoring the rest deconfounds correlated
#' regions: a region masked by a stronger, correlated neighbour resurfaces once
#' that neighbour has been taken out. This typically gives much better recall
#' than a single one-shot ranking, at a fraction of the cost (no model fit).
#'
#' @param y numeric response vector of length `n`.
#' @param G_list length-`p` list of `n x n` per-region covariance (GRM)
#'   matrices. Ignored when `V_list` is supplied.
#' @param K how many leading eigenvectors to keep per region. `"mp"` (default)
#'   chooses each region's count adaptively from a Marchenko-Pastur eigenvalue
#'   threshold (keep eigenvalues above `median(ev) * (1 + sqrt(r/n))^2`); a
#'   fraction in `(0,1)` keeps the top directions explaining that share of the
#'   region's variance; a positive integer fixes the count. Capped per region
#'   at its number of positive eigenvalues.
#' @param max_select optional hard cap on how many regions to select. `NULL`
#'   (default) means no artificial cap: selection is driven entirely by the
#'   noise-floor stop (`stop_rule = "noise"`), so it adapts to however many
#'   regions actually carry signal. Set an integer only to bound the refit
#'   budget (or when using `stop_rule = "fixed"`, where it sets the count).
#' @param standardize if TRUE (default) centre and scale `y` to unit variance
#'   before selection. Recommended: it makes the noise-floor stopping rule
#'   scale-free. (Region scaling is irrelevant to OMP, which uses eigenvectors
#'   only, so kernels do not need standardising.)
#' @param stop_rule stopping rule. `"noise"` (default) stops early when the
#'   best remaining score falls to the residual noise floor; `"fixed"` always
#'   selects `max_select` regions.
#' @param noise_mult multiplier on the noise floor for the `"noise"` rule
#'   (default 1.5, recall-first; larger = stop sooner, fewer regions selected).
#' @param union_marginal if TRUE (default) add, to the greedy selection, any
#'   region whose *marginal* score (projection onto the original response, no
#'   residualisation) clears the noise floor --- a safety net that recovers a
#'   strong region the greedy step may have "removed away" via a correlated
#'   neighbour. Capped so the candidate set does not balloon; the refit shrinks
#'   the harmless extras.
#' @param V_list optional precomputed list of `n x K` orthonormal eigenvector
#'   matrices (one per region). When supplied, `G_list` is not needed and the
#'   (expensive) eigendecompositions are skipped --- useful for repeated calls
#'   on the same set of regions (e.g. many phenotypes / replicates).
#' @param verbose logical; if TRUE, report each pick.
#' @return An object of class `"omp_select"`: a list with
#'   \describe{
#'     \item{selected}{region indices in selection order (best first).}
#'     \item{scores}{the score at the moment each region was selected.}
#'     \item{n_selected}{number of regions selected.}
#'     \item{residual_var}{`||r||^2 / n` after each pick (the shrinking residual).}
#'     \item{V_list}{the per-region eigenvectors, so a later call can reuse them
#'       via the `V_list` argument.}
#'   }
#' @seealso [cblocks()], [pblocks()] for the exact refit on the selected regions.
#' @examples
#' \dontrun{
#' ## Stage 1: cheap OMP selection
#' sel <- omp_select(y, G_list, max_select = 10)
#' k   <- sel$selected                      # candidate regions, best first
#'
#' ## Stage 2: exact refit with cblocks on just those regions
#' cm  <- cblocks(length(y), length(k), G_list[k], prior = "horseshoe")
#' fit <- INLA::inla(y ~ 1 + f(idx, model = cm),
#'                   data = data.frame(y = y, idx = seq_along(y)),
#'                   family = "gaussian")
#' ## per-region variance estimates (back on the y scale):
#' v <- exp(-fit$summary.hyperpar$mode[-1]) * var(y)
#' }
#' @export
omp_select <- function(y, G_list = NULL, K = "mp", max_select = NULL,
                       standardize = TRUE, stop_rule = c("noise", "fixed"),
                       noise_mult = 1.5, union_marginal = TRUE,
                       V_list = NULL, verbose = FALSE) {
  stop_rule <- match.arg(stop_rule)
  if (is.null(max_select) && stop_rule == "fixed")
    stop('`stop_rule = "fixed"` requires an explicit `max_select`', call. = FALSE)
  y <- as.numeric(y); n <- length(y)
  if (is.character(K) && !identical(K, "mp"))
    stop('character `K` must be "mp"', call. = FALSE)
  if (standardize) {
    y <- y - mean(y)
    s <- sqrt(sum(y^2) / n)
    if (s > 0) y <- y / s
  }

  ## --- build / validate per-region orthonormal eigenbases ---
  if (is.null(V_list)) {
    if (is.null(G_list))
      stop("supply either `G_list` or `V_list`", call. = FALSE)
    p <- length(G_list)
    V_list <- vector("list", p)
    for (j in seq_len(p)) {
      Gj <- as.matrix(G_list[[j]])
      if (nrow(Gj) != n || ncol(Gj) != n)
        stop(sprintf("G_list[[%d]] must be %d x %d", j, n, n), call. = FALSE)
      e  <- eigen((Gj + t(Gj)) / 2, symmetric = TRUE)
      Kj <- .omp_choose_K(e$values, n, K)                # adaptive per-region rank
      V_list[[j]] <- e$vectors[, seq_len(Kj), drop = FALSE]
    }
  } else {
    p <- length(V_list)
    for (j in seq_len(p))
      if (nrow(V_list[[j]]) != n)
        stop(sprintf("V_list[[%d]] must have %d rows", j, n), call. = FALSE)
  }

  ## --- OMP loop ---
  max_select <- if (is.null(max_select)) p else min(as.integer(max_select), p)
  r <- y
  remaining <- seq_len(p)
  selected <- integer(0); scores <- numeric(0); resid_var <- numeric(0)
  for (step in seq_len(max_select)) {
    sc <- vapply(remaining, function(j) {
      Vj <- V_list[[j]]
      sum(crossprod(Vj, r)^2) / ncol(Vj)
    }, numeric(1))
    bi <- which.max(sc)
    best <- remaining[bi]; best_sc <- sc[bi]
    if (stop_rule == "noise" && best_sc < noise_mult * sum(r^2) / n) break
    Vb <- V_list[[best]]
    r <- r - Vb %*% crossprod(Vb, r)            # remove region's contribution
    selected   <- c(selected, best)
    scores     <- c(scores, best_sc)
    resid_var  <- c(resid_var, sum(r^2) / n)
    remaining  <- remaining[-bi]
    if (verbose)
      message(sprintf("step %2d: region %d  (score %.4g, residual var %.4g)",
                      step, best, best_sc, sum(r^2) / n))
  }

  ## --- union with marginal ranking (recall safety net) ---
  ## The greedy loop can "remove away" a strong region that is correlated with a
  ## selected neighbour. But such a region keeps a high MARGINAL score (projection
  ## onto the ORIGINAL response, no residualisation), so it cannot be hidden. We
  ## add any region whose marginal score clears the noise floor, capped so the
  ## candidate set does not balloon (the cblocks refit shrinks the harmless extras).
  if (isTRUE(union_marginal)) {
    marg <- vapply(seq_len(p), function(j) {
      Vj <- V_list[[j]]; sum(crossprod(Vj, y)^2) / ncol(Vj)
    }, numeric(1))
    floor0 <- noise_mult * sum(y^2) / n
    extra  <- setdiff(which(marg > floor0), selected)
    extra  <- extra[order(marg[extra], decreasing = TRUE)]
    n_omp  <- length(selected)
    cap    <- min(max_select, n_omp + max(5L, n_omp))     # soft cap: don't balloon
    room   <- cap - n_omp
    if (room > 0 && length(extra) > 0) {
      add      <- extra[seq_len(min(room, length(extra)))]
      selected <- c(selected, add)
      scores   <- c(scores, marg[add])
    }
  }

  structure(list(selected = selected, scores = scores,
                 n_selected = length(selected), residual_var = resid_var,
                 V_list = V_list, call = match.call()),
            class = "omp_select")
}

## Choose how many leading eigenvectors to keep for one region, given its
## eigenvalues `ev` (descending), sample size `n`, and the rule `K`:
##   "mp"          -> Marchenko-Pastur edge: keep eigenvalues above
##                    median(ev) * (1 + sqrt(r/n))^2   (r = # positive eigenvalues)
##   fraction (0,1)-> keep top directions explaining that share of sum(ev)
##   integer >= 1  -> fixed count
## Ported from the former rblocks selector (.rblocks_choose_K).
.omp_choose_K <- function(ev, n, K) {
  ev <- ev[ev > ev[1] * 1e-8]               # positive (non-negligible) eigenvalues
  r  <- length(ev)
  if (r == 0L) return(0L)
  if (is.character(K)) {                    # "mp"
    edge <- stats::median(ev) * (1 + sqrt(r / n))^2
    return(max(1L, sum(ev > edge)))
  }
  if (K > 0 && K < 1) {                     # variance fraction (inertia)
    cum <- cumsum(ev) / sum(ev)
    return(max(1L, which(cum >= K)[1]))
  }
  min(as.integer(K), r)                     # fixed count
}

#' @export
print.omp_select <- function(x, ...) {
  cat(sprintf("OMP selection: %d region(s)\n", x$n_selected))
  if (x$n_selected > 0L) {
    cat("  order :", paste(x$selected, collapse = " "), "\n")
    cat("  scores:", paste(sprintf("%.3g", x$scores), collapse = " "), "\n")
  }
  invisible(x)
}
