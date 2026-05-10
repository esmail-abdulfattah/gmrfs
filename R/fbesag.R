## Partitioned Besag model.

.fbesag_pack <- function(graph, id_s, sd_gamma, param) {
  lam <- -log(param$p2) / param$p1
  n   <- dim(graph)[1]
  P   <- length(unique(id_s))
  g   <- INLA::inla.read.graph(graph)

  constr_inter <- list(A = matrix(1, 1, dim(graph)[1]), e = rep(0, 1))
  scaled_graph <- as.matrix(INLA:::inla.scale.model(graph, constr_inter))
  scaled_cnst  <- scaled_graph[1, 1] / graph[1, 1]

  sigm2 <- sd_gamma * sd_gamma
  e <- eigen(diag(P) - (1 / P) * matrix(1, P, P))
  D <- diag(c(1.0 / e$values[1:(P - 1)]))
  inv_tilda_Sigma <- (1 / sigm2) * e$vectors[, 1:(P - 1)] %*% D %*% t(e$vectors[, 1:(P - 1)])
  val1 <- inv_tilda_Sigma[1, 1]
  val2 <- inv_tilda_Sigma[1, 2]

  graph_vec <- local({
    ii <- integer(0); jj <- integer(0)
    for (i in 1:g$n) {
      ind <- which(g$nbs[[i]] >= i)
      if (length(ind) > 0) {
        new_n <- g$nbs[[i]][ind]
        ii <- c(ii, i, rep(i, length(new_n)))
        jj <- c(jj, i, new_n)
      } else {
        ii <- c(ii, i); jj <- c(jj, i)
      }
    }
    c(g$n, length(ii), ii - 1, jj - 1)
  })
  graph_vec <- c(length(graph_vec), graph_vec)

  misc_vec <- local({
    one_vector <- integer(0)
    for (i in 1:g$n) {
      num_nei_i <- g$nnbs[i]
      one_vector <- c(one_vector, num_nei_i, id_s[i] - 1)
      g$nbs[[i]] <- sort(stats::na.omit(g$nbs[[i]]))
      size_neighbors <- length(g$nbs[[i]])
      one_vector <- c(one_vector, size_neighbors)
      for (j in seq_len(size_neighbors)) {
        tick <- g$nbs[[i]][j]
        one_vector <- c(one_vector, id_s[tick] - 1, tick - 1)
      }
    }
    one_vector
  })

  list(graph_vec = graph_vec,
       misc_vec  = misc_vec,
       lam       = lam,
       P         = P,
       n         = n,
       invSig    = c(val1, val2, scaled_cnst))
}

#' Partitioned Besag model (fbesag)
#'
#' @param graph adjacency matrix or anything `INLA::inla.read.graph` accepts.
#' @param id integer vector (length = nrow(graph)) of partition ids in 1..P.
#' @param sd_gamma prior sd for partition effects (default 0.2).
#' @param param list with `p1`, `p2` defining lambda = -log(p2)/p1.
#' @param initial length-P vector of starting log-precisions, or `-999`
#'   sentinel for the default (rep(4, P)).
#' @param shlib optional path to a specific `.so`/`.dll`. When supplied it
#'   overrides the resolution chain (local dir, INLA, package).
#' @param prefer_inla if TRUE (default) and INLA ships `libfbesag.so`, use it;
#'   otherwise fall back to this package's compiled library.
#' @return cgeneric model handle for use as `f(idx, model = .)` in INLA.
#' @export
fbesag <- function(graph, id,
                         sd_gamma    = 0.2,
                         param       = list(p1 = 1, p2 = 1e-5),
                         initial     = -999,
                         shlib       = NULL,
                         prefer_inla = TRUE) {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("INLA is required: install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable')",
         call. = FALSE)

  res <- .fbesag_pack(graph, id, sd_gamma, param)

  if (length(initial) == 1L && initial[1] == -999) {
    initial <- rep(4, res$P)
  } else if (length(initial) != res$P) {
    stop(sprintf("initial must have length P = %d", res$P), call. = FALSE)
  }

  INLA::inla.cgeneric.define(
    model              = "inla_cgeneric_fbesag",
    shlib              = gmrfs_resolve_shlib("fbesag", shlib, prefer_inla),
    n                  = as.integer(res$n),
    npart              = as.integer(res$P),
    VEC_CGENERIC_GRAPH = as.integer(res$graph_vec),
    debug              = FALSE,
    lam                = c(res$lam),
    invSig             = res$invSig,
    misc               = as.integer(res$misc_vec),
    initial            = as.numeric(initial)
  )
}
