## Partitioned Besag model.

.fbesag_pack <- function(graph, id_s, sd_gamma, param) {
  lam <- -log(param$p2) / param$p1
  P   <- length(unique(id_s))
  if (FALSE) g <- INLA::inla.read.graph(graph)
  g   <- gmrfs_read_graph(graph)
  n   <- g$n
  if (length(id_s) != n || !all(sort(unique(id_s)) == seq_len(P)))
    stop("`id` must have one entry per node, with values 1..P all used", call. = FALSE)

  ## inla.read.graph keeps the order of a graph file; the graph and the Q
  ## values below must list the neighbours in the same order, so sort once
  g$nbs <- lapply(g$nbs, function(nb) sort(as.integer(nb)))

  if (FALSE) {
    ## old: one scaling constant from node 1, assumes a connected graph
    ## and that `graph` is the structure matrix R = D - W
    constr_inter <- list(A = matrix(1, 1, dim(graph)[1]), e = rep(0, 1))
    scaled_graph <- as.matrix(INLA::inla.scale.model(graph, constr_inter))
    scaled_cnst  <- scaled_graph[1, 1] / graph[1, 1]
  }

  ## new: as inla_besag_scale() in INLA. One scaling constant per connected
  ## component (geometric mean of the marginal variances of R_k under its own
  ## sum-to-zero constraint). Islands get the code -1 and are treated
  ## specially in the C code (Q_ii = tau_p(i)).
  if (FALSE) R <- -INLA::inla.graph2matrix(g)
  R <- -gmrfs_graph2matrix(g)
  diag(R) <- 0
  diag(R) <- -Matrix::rowSums(R)
  cc_nodes <- g$cc$nodes
  cc_size  <- lengths(cc_nodes)
  scale_vec <- numeric(n)
  for (k in seq_along(cc_nodes)) {
    idx <- cc_nodes[[k]]
    if (cc_size[k] == 1L) {
      scale_vec[idx] <- -1
    } else {
      Rk <- R[idx, idx, drop = FALSE]
      if (FALSE) Sk <- INLA::inla.scale.model(Rk, list(A = matrix(1, 1, length(idx)), e = 0))
      Sk <- gmrfs_scale_model(Rk, list(A = matrix(1, 1, length(idx)), e = 0))
      scale_vec[idx] <- Sk[1, 1] / Rk[1, 1]
    }
  }
  scaled_cnst <- scale_vec[1]   # kept for the old C code path

  ## one sum-to-zero constraint per connected component of size >= 2
  cc_big <- which(cc_size >= 2L)
  A_cc <- matrix(0, length(cc_big), n)
  for (k in seq_along(cc_big)) A_cc[k, cc_nodes[[cc_big[k]]]] <- 1

  ## besag's normalising constant leaves out 0.5 * log pdet of the scaled
  ## structure matrix; INLA's own computation for this cgeneric keeps it. K
  ## removes it, so with one partition fbesag has besag's mlik. Per component,
  ## pdet(c_k R_k) = c_k^(n_k - 1) * n_k * det(R_k[-1, -1]) (matrix-tree);
  ## islands contribute log 1 = 0.
  log_pdet <- 0
  for (k in cc_big) {
    idx <- cc_nodes[[k]]
    Rk <- R[idx, idx, drop = FALSE]
    log_pdet <- log_pdet + (length(idx) - 1) * log(scale_vec[idx[1]]) + log(length(idx)) +
      as.numeric(Matrix::determinant(Rk[-1, -1, drop = FALSE], logarithm = TRUE)$modulus)
  }
  K <- -0.5 * log_pdet

  sigm2 <- sd_gamma * sd_gamma
  if (FALSE) {
    ## old: fails for P = 1 (1:(P - 1) is 1:0)
    e <- eigen(diag(P) - (1 / P) * matrix(1, P, P))
    D <- diag(c(1.0 / e$values[1:(P - 1)]))
    inv_tilda_Sigma <- (1 / sigm2) * e$vectors[, 1:(P - 1)] %*% D %*% t(e$vectors[, 1:(P - 1)])
    val1 <- inv_tilda_Sigma[1, 1]
    val2 <- inv_tilda_Sigma[1, 2]
  }
  ## pseudo-inverse of sigm2 (I - J/P) is (I - J/P) / sigm2; no deviations if P = 1
  val1 <- if (P > 1) (1 - 1 / P) / sigm2 else 0
  val2 <- if (P > 1) -(1 / P) / sigm2 else 0

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
       invSig    = c(val1, val2, scaled_cnst),
       scale_vec = scale_vec,
       A_cc      = A_cc,
       K         = K)
}

#' Partitioned Besag model (fbesag)
#'
#' @param graph adjacency matrix, structure matrix, `inla.graph`, spdep `nb`
#'   object or graph file (what `INLA::inla.read.graph()` accepts). INLA is
#'   not needed to build the model.
#' @param id integer vector (length = nrow(graph)) of partition ids in 1..P.
#' @param sd_gamma prior sd for partition effects (default 0.2).
#' @param param list with `p1`, `p2` defining lambda = -log(p2)/p1.
#' @param initial length-P vector of starting log-precisions, or `-999`
#'   sentinel for the default (rep(4, P)).
#' @param shlib optional path to a specific `.so`/`.dll`. When supplied it
#'   overrides the resolution chain (local dir, INLA, package).
#' @param prefer_inla if TRUE (default) and INLA ships `libfbesag.so`, use it;
#'   otherwise fall back to this package's compiled library.
#' @param constr if TRUE (default), add one sum-to-zero constraint per
#'   connected component of size >= 2 (islands are proper and get none).
#' @param extraconstr optional `list(A, e)` of further constraints; appended
#'   to the per-component ones.
#' @param nrep number of replicates when used with `f(..., replicate = )`.
#'   The constant that puts the mlik on besag's convention is per replicate,
#'   and the model cannot see `replicate`, so give it here (default 1).
#' @return cgeneric model handle for use as `f(idx, model = .)` in INLA. The
#'   constraints and `rankdef` are carried in the handle, so `f()` picks them
#'   up unless overridden there.
#' @examples
#' ## path graph on 5 nodes, two partitions
#' W <- matrix(0, 5, 5)
#' W[cbind(1:4, 2:5)] <- 1
#' W <- W + t(W)
#' fm <- fbesag(W, id = c(1, 1, 1, 2, 2))
#' gmrf_Q(fm)
#' @export
fbesag <- function(graph, id,
                         sd_gamma    = 0.2,
                         param       = list(p1 = 1, p2 = 1e-5),
                         initial     = -999,
                         shlib       = NULL,
                         prefer_inla = TRUE,
                         constr      = TRUE,
                         extraconstr = NULL,
                         nrep        = 1L) {
  if (FALSE) {
    if (!requireNamespace("INLA", quietly = TRUE))
      stop("INLA is required: install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable')",
           call. = FALSE)
  }

  res <- .fbesag_pack(graph, id, sd_gamma, param)

  if (length(initial) == 1L && initial[1] == -999) {
    initial <- rep(4, res$P)
  } else if (length(initial) != res$P) {
    stop(sprintf("initial must have length P = %d", res$P), call. = FALSE)
  }

  ## INLA::inla.cgeneric.define() if INLA is installed, else the same object
  ## built here
  cm <- gmrfs_cgeneric_define(
    "inla_cgeneric_gmrfs_fbesag", "fbesag",
    shlib              = shlib,
    prefer_inla        = prefer_inla,
    n                  = as.integer(res$n),
    npart              = as.integer(res$P),
    VEC_CGENERIC_GRAPH = as.integer(res$graph_vec),
    debug              = FALSE,
    lam                = c(res$lam),
    invSig             = res$invSig,
    misc               = as.integer(res$misc_vec),
    initial            = as.numeric(initial),
    scale              = as.numeric(res$scale_vec),
    sd_gamma           = as.numeric(sd_gamma),
    K                  = as.numeric(nrep * res$K)
  )

  ## constraints, as f() does for besag/bym2 with adjust.for.con.comp = TRUE.
  ## f() copies model$f$<arg> into its own argument when that is missing.
  A <- if (isTRUE(constr)) res$A_cc else matrix(0, 0, res$n)
  e <- rep(0, nrow(A))
  if (!is.null(extraconstr)) {
    extraconstr$A <- matrix(extraconstr$A, ncol = res$n)
    A <- rbind(A, extraconstr$A)
    e <- c(e, extraconstr$e)
  }
  if (nrow(A) > 0) cm$f$extraconstr <- list(A = A, e = e)
  ## as in f(): rankdef = #components of size >= 2 (+ rows of extraconstr)
  cm$f$rankdef <- nrow(res$A_cc) + (if (is.null(extraconstr)) 0L else nrow(extraconstr$A))
  ## INLA adds diagonal = 1e-4 to Q of a model with constraints. INLA computes
  ## this model's normalising constant from Q, so the 1e-4 enters the mlik
  ## (0.26 too high on a 150 x 150 lattice); with 1e-8 the mlik is exact to
  ## 1e-3 on that lattice. A diagonal given in f() still takes precedence.
  cm$f$diagonal <- 1e-8
  cm
}
