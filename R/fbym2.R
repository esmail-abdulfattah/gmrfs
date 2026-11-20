## Flexible BYM2: bym2 whose structured part is fbesag.

.fbym2_pack <- function(graph, id, phi) {
  if (FALSE) g <- INLA::inla.read.graph(graph)
  g <- gmrfs_read_graph(graph)
  n <- g$n
  P <- length(unique(id))
  if (length(id) != n || !all(sort(unique(id)) == seq_len(P)))
    stop("`id` must have one entry per node, with values 1..P all used", call. = FALSE)

  ## structure matrix R = D - W and one scaling constant per connected
  ## component (as inla_besag_scale); islands get the code -1
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

  ## sorted neighbour lists, used by both the graph and misc vectors
  nbs <- lapply(seq_len(n), function(i) sort(g$nbs[[i]]))

  ## graph of x = (v, u), upper triangle, in the order the C code fills Q:
  ## rows 1..n: (i, i), (i, i+n); rows n+1..2n: (n+i, n+i), (n+i, n+j), j > i
  up <- lapply(seq_len(n), function(i) nbs[[i]][nbs[[i]] > i])
  ii <- c(rep(seq_len(n), each = 2),
          unlist(lapply(seq_len(n), function(i) rep(n + i, 1L + length(up[[i]])))))
  jj <- c(as.vector(rbind(seq_len(n), n + seq_len(n))),
          unlist(lapply(seq_len(n), function(i) c(n + i, n + up[[i]]))))
  graph_vec <- c(2L * n, length(ii), ii - 1L, jj - 1L)
  graph_vec <- c(length(graph_vec), graph_vec)

  ## per node: #neighbours, partition, #neighbours, then (partition, index) pairs
  misc_vec <- unlist(lapply(seq_len(n), function(i) {
    nb <- nbs[[i]]
    c(length(nb), id[i] - 1L, length(nb), as.vector(rbind(id[nb] - 1L, nb - 1L)))
  }))

  ## bym2 PC prior on phi, tabulated on the logit scale
  tab <- gmrfs_pc_bym_phi(graph = g, u = phi[1], alpha = phi[2], return.as.table = TRUE)
  tail <- attr(tab, "tail")
  tab <- as.numeric(strsplit(trimws(sub("^table:", "", tab)), "[ ]+")[[1]])
  m <- length(tab) / 2
  phi_x <- tab[seq_len(m)]
  phi_y <- tab[m + seq_len(m)]

  ## Beyond the table (logit phi > theta_max) the prior has the exact tail of
  ## the PC prior: d(theta)^2 = d_max^2 + r (L(theta) - L(theta_max)), with
  ## L(theta) = log(1 + e^theta), and log prior c0 - lambda d + log(r phi / d).
  ## c0 makes it continuous with the table. Inside the table nothing changes,
  ## so fbym2 keeps bym2's prior there; INLA's bym2 keeps the prior constant
  ## beyond the table, which is not a proper density.
  phi_tail <- c(tail[["theta_max"]], tail[["d_max"]], tail[["rankdef"]], tail[["lambda"]], 0)
  phi_tail[5] <- phi_y[m] + tail[["lambda"]] * tail[["d_max"]] -
    log(tail[["rankdef"]] * stats::plogis(tail[["theta_max"]]) / tail[["d_max"]])
  ## total mass of the prior on the whole line, for the exact (normalised) prior
  ## in gmrf_log_prior(): the table (trapezoid), the linear left tail and the
  ## right tail, int exp(c0 - lambda d) 2 dd = 2 exp(c0 - lambda d_max) / lambda
  slope_left <- (phi_y[2] - phi_y[1]) / (phi_x[2] - phi_x[1])
  mass <- sum(diff(phi_x) * (exp(phi_y[-1]) + exp(phi_y[-m])) / 2) +
    exp(phi_y[1]) / slope_left +
    2 * exp(phi_tail[5] - tail[["lambda"]] * tail[["d_max"]]) / tail[["lambda"]]

  ## bym2's normalising constant leaves out 0.5 * log pdet(R*) and the
  ## log(2 pi) terms of the n - r constrained spatial dimensions; INLA's own
  ## computation for this cgeneric keeps them. K removes them, so fbym2 has
  ## bym2's mlik (and with one partition the same value). Per component,
  ## pdet(c_k R_k) = c_k^(n_k - 1) * n_k * det(R_k[-1, -1]) (matrix-tree);
  ## islands contribute log 1 = 0.
  cc_big <- which(cc_size >= 2L)
  log_pdet <- 0
  for (k in cc_big) {
    idx <- cc_nodes[[k]]
    Rk <- R[idx, idx, drop = FALSE]
    log_pdet <- log_pdet + (length(idx) - 1) * log(scale_vec[idx[1]]) + log(length(idx)) +
      as.numeric(Matrix::determinant(Rk[-1, -1, drop = FALSE], logarithm = TRUE)$modulus)
  }
  K <- -0.5 * log_pdet + 0.5 * (n - length(cc_big)) * log(2 * pi)

  ## one sum-to-zero constraint per component of size >= 2, on the u half
  A_cc <- matrix(0, length(cc_big), 2L * n)
  for (k in seq_along(cc_big)) A_cc[k, n + cc_nodes[[cc_big[k]]]] <- 1

  list(n = n, P = P, graph_vec = graph_vec, misc_vec = misc_vec,
       scale_vec = scale_vec, phi_x = phi_x, phi_y = phi_y, phi_tail = phi_tail,
       phi_log_mass = log(mass), A_cc = A_cc, K = K)
}

#' Flexible BYM2 model (fbym2)
#'
#' BYM2 whose structured part is the partitioned Besag (fbesag):
#' `v = (sqrt(1 - phi) w + sqrt(phi) u) / sqrt(tau)`, with `w` iid and `u` an
#' fbesag field with relative partition precisions `exp(d_k)`, `sum(d) = 0`.
#' The iid part sees only the global `tau`; the spatial part of partition `k`
#' has precision `tau * exp(d_k)`. With all `d = 0` this is INLA's `bym2`.
#'
#' The latent vector has length `2n`: entries `1..n` are `v` (the total effect,
#' used in the predictor), entries `n+1..2n` are the standardised spatial
#' part `u`. Hyperparameters: `log tau`, `logit phi`, `d_1, ..., d_(P-1)`.
#'
#' The marginal likelihood follows INLA's `bym2` convention (the graph-only
#' terms `bym2` leaves out of its normalising constant are left out here too),
#' so `mlik` can be compared with `bym2` directly; with one partition it is
#' the same.
#'
#' @param graph adjacency matrix, structure matrix, `inla.graph`, spdep `nb`
#'   object or graph file. INLA is not needed to build the model.
#' @param id integer vector (one per node) of partition ids in 1..P.
#' @param sd_gamma prior sd of the partition deviations `d` (default 0.2).
#' @param prec PC prior `c(u, alpha)` on the precision `tau`,
#'   `P(1/sqrt(tau) > u) = alpha` (default as bym2: `c(1, 0.01)`).
#' @param phi PC prior `c(u, alpha)` on `phi`, `P(phi < u) = alpha`
#'   (default as bym2: `c(0.5, 0.5)`). The prior is tabulated for `logit(phi)`
#'   up to 12, with the same values as INLA's `bym2`; beyond, it has the
#'   exact tail of the PC prior (INLA's `bym2` keeps it constant there).
#' @param initial initial values of the hyperparameters (length `2 + P - 1`);
#'   default `c(4, -3, 0, ...)`.
#' @param shlib optional path to a specific `.so`/`.dll`.
#' @param prefer_inla if TRUE (default) and INLA ships `libfbym2.so`, use it.
#' @param constr if TRUE (default), one sum-to-zero constraint on `u` per
#'   connected component of size >= 2.
#' @param extraconstr optional `list(A, e)` of further constraints on the
#'   `2n` latent vector; appended to the per-component ones.
#' @param nrep number of replicates when used with `f(..., replicate = )`.
#'   The constant that puts the mlik on bym2's convention is per replicate,
#'   and the model cannot see `replicate`, so give it here (default 1).
#' @return cgeneric model handle for use as `f(idx, model = .)` in INLA, with
#'   the constraints and `rankdef` carried in the handle.
#' @examples
#' \donttest{
#' ## path graph on 5 nodes, two partitions
#' W <- matrix(0, 5, 5)
#' W[cbind(1:4, 2:5)] <- 1
#' W <- W + t(W)
#' op <- options(mc.cores = 1)   # the prior of phi is tabulated in parallel
#' bm <- fbym2(W, id = c(1, 1, 1, 2, 2))
#' options(op)
#' dim(gmrf_Q(bm))
#' }
#' @export
fbym2 <- function(graph, id,
                  sd_gamma    = 0.2,
                  prec        = c(1, 0.01),
                  phi         = c(0.5, 0.5),
                  initial     = NULL,
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

  res <- .fbym2_pack(graph, id, phi)
  ntheta <- 2L + res$P - 1L
  if (is.null(initial)) {
    initial <- c(4, -3, rep(0, res$P - 1L))
  } else if (length(initial) != ntheta) {
    stop(sprintf("initial must have length 2 + P - 1 = %d", ntheta), call. = FALSE)
  }

  ## INLA::inla.cgeneric.define() if INLA is installed, else the same object
  ## built here
  cm <- gmrfs_cgeneric_define(
    "inla_cgeneric_gmrfs_fbym2", "fbym2",
    n        = as.integer(2L * res$n),
    npart    = as.integer(res$P),
    graph    = as.integer(res$graph_vec),
    debug    = FALSE,
    shlib    = shlib,
    prefer_inla = prefer_inla,
    misc     = as.integer(res$misc_vec),
    scale    = as.numeric(res$scale_vec),
    prec     = as.numeric(prec),
    phi_x    = as.numeric(res$phi_x),
    phi_y    = as.numeric(res$phi_y),
    phi_tail = as.numeric(res$phi_tail),
    phi_log_mass = as.numeric(res$phi_log_mass),
    sd_gamma = as.numeric(sd_gamma),
    initial  = as.numeric(initial),
    K        = as.numeric(nrep * res$K)
  )

  ## constraints, as f() does for bym2: on the u half, one per component
  A <- if (isTRUE(constr)) res$A_cc else matrix(0, 0, 2L * res$n)
  e <- rep(0, nrow(A))
  if (!is.null(extraconstr)) {
    extraconstr$A <- matrix(extraconstr$A, ncol = 2L * res$n)
    A <- rbind(A, extraconstr$A)
    e <- c(e, extraconstr$e)
  }
  if (nrow(A) > 0) cm$f$extraconstr <- list(A = A, e = e)
  cm$f$rankdef <- nrow(res$A_cc) + (if (is.null(extraconstr)) 0L else nrow(extraconstr$A))
  ## INLA adds diagonal = 1e-4 to Q of a model with constraints. INLA computes
  ## this model's normalising constant from Q, so the 1e-4 enters the mlik
  ## (0.26 too high on a 150 x 150 lattice); with 1e-8 the mlik is exact to
  ## 1e-3 on that lattice. A diagonal given in f() still takes precedence.
  cm$f$diagonal <- 1e-8
  cm
}
