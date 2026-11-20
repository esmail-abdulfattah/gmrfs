## The PC prior for the mixing parameter phi of the bym2 model.
##
## Copied from R-INLA (rinla/R/pc-bym.R, commit f711ebb82), where
## inla.pc.bym.phi() is not exported; renamed. INLA's helpers are replaced by
## the copies in graph.R, and the constrained log-determinant, which INLA gets
## from its binary (inla.qsample()), is computed from a sparse Cholesky
## factorisation, so INLA is not needed. R-INLA is Copyright Haavard Rue,
## MIT licence.

## inla.mclapply() without the INLA options: `mc.cores` if set, else all
## physical cores
gmrfs_mclapply <- function(...) {
  if (.Platform$OS.type == "windows") return(lapply(...))
  cores <- getOption("mc.cores", parallel::detectCores(logical = FALSE))
  if (is.na(cores)) cores <- 1L
  parallel::mclapply(..., mc.cores = cores)
}

gmrfs_bym_constr <- function(Q, adjust.for.con.comp = TRUE) {
  ## return the rankdef and the constr for the BYM model given in Q

  n <- dim(Q)[1]
  if (FALSE) g <- INLA::inla.read.graph(Q)
  g <- gmrfs_read_graph(Q)
  cc.n <- lengths(g$cc$nodes)
  cc.n1 <- sum(cc.n == 1L)
  cc.n2 <- sum(cc.n >= 2L)

  if (adjust.for.con.comp) {
    constr <- list(A = matrix(0, cc.n2, n), e = rep(0, cc.n2))
    kk <- 0
    for (k in which(cc.n >= 2L)) {
      kk <- kk + 1
      constr$A[kk, g$cc$nodes[[k]]] <- 1
    }
    stopifnot(kk == cc.n2)
    rankdef <- cc.n2
  } else {
    nc <- 1
    constr <- list(A = matrix(1, nc, n), e = rep(1, nc))
    rankdef <- cc.n1 + 1
  }

  return(list(
    constr = constr,
    rankdef = rankdef,
    cc.n = cc.n,
    cc.n1 = cc.n1,
    cc.n2 = cc.n2
  ))
}

gmrfs_sparse_det_bym <- function(Q,
                rankdef,
                adjust.for.con.comp = TRUE,
                log = TRUE,
                constr = NULL,
                eps = 0.01 * sqrt(.Machine$double.eps)) {
  ## compute the (log-)determinant. Assume a sum to zero constraint on each connected
  ## component >1. if 'constr' is given,  it is supposed to be the 'constr' that is
  ## constructed below.

  adjust.for.con.comp <- as.logical(adjust.for.con.comp)
  stopifnot(adjust.for.con.comp)
  if (FALSE) Q <- INLA::inla.as.sparse(Q)
  Q <- gmrfs_as_sparse(Q)
  n <- dim(Q)[1]
  res <- NULL

  if (is.null(constr)) {
    res <- gmrfs_bym_constr(Q, adjust.for.con.comp = adjust.for.con.comp)
    constr <- res$constr
  }
  if (missing(rankdef)) {
    if (!is.null(constr)) {
      rankdef <- nrow(constr$A)
    } else {
      rankdef <- 0
    }
  }
  d <- diag(Q)
  diag(Q) <- d + mean(d) * eps
  if (FALSE) {
    res <- INLA::inla.qsample(n = 1, Q = Q, sample = matrix(0, n, 1), constr = constr)
    logdet <- 2.0 * (res$logdens + (n - rankdef) / 2.0 * log(2.0 * pi))
  }
  ## twice the log density at x = 0 of N(0, Q^-1) given A x = 0, plus
  ## (n - rankdef) log(2 pi): log|Q| + log|A Q^-1 A'| - log|A A'|. The last
  ## term does not depend on Q and cancels in the prior of phi.
  ch <- Matrix::Cholesky(Matrix::forceSymmetric(methods::as(Q, "CsparseMatrix"), uplo = "U"),
                         LDL = FALSE, perm = TRUE, super = FALSE)
  logdet <- as.numeric(Matrix::determinant(ch, logarithm = TRUE, sqrt = FALSE)$modulus)
  if (!is.null(constr) && NROW(constr$A) > 0) {
    A <- matrix(constr$A, ncol = n)
    logdet <- logdet + as.numeric(determinant(A %*% as.matrix(Matrix::solve(ch, t(A))))$modulus) -
      as.numeric(determinant(tcrossprod(A))$modulus)
  }
  return(if (log) logdet else exp(logdet))
}

gmrfs_pc_bym_Q <- function(graph) {
  if (FALSE) Q <- -INLA::inla.graph2matrix(graph)
  Q <- -gmrfs_graph2matrix(graph)
  diag(Q) <- 0
  diag(Q) <- -Matrix::rowSums(Q)
  return(Q)
}

gmrfs_pc_bym_phi <- function(graph,
              Q,
              eigenvalues = NULL,
              marginal.variances = NULL,
              rankdef,
              alpha,
              u = 1 / 2,
              lambda,
              scale.model = TRUE,
              return.as.table = FALSE,
              adjust.for.con.comp = TRUE,
              ## where to switch to alternative strategy
              eps = sqrt(.Machine$double.eps),
              debug = FALSE) {
  my.debug <- function(...) if (debug) cat("*** debug *** gmrfs_pc_bym_phi: ", ..., "\n")

  scale.model <- as.logical(scale.model)
  adjust.for.con.comp <- as.logical(adjust.for.con.comp)

  ## I must assume this!!!
  stopifnot(scale.model)
  stopifnot(adjust.for.con.comp)
  res <- NULL

  if (missing(eigenvalues) && missing(marginal.variances)) {
    ## computes the prior for the mixing parameter `phi'.
    if (missing(graph) && !missing(Q)) {
      Q <- gmrfs_as_sparse(Q)
    } else if (!missing(graph) && missing(Q)) {
      Q <- gmrfs_pc_bym_Q(graph)
    } else if (missing(graph) && missing(Q)) {
      stop("Either <Q> or <graph> must be given.")
    } else {
      stop("Only one of <Q> and <graph> can be given.")
    }

    res <- gmrfs_bym_constr(Q, adjust.for.con.comp = adjust.for.con.comp)
    if (missing(rankdef)) {
      rankdef <- res$rankdef
    }

    n <- dim(Q)[1]
    if (FALSE) res <- INLA::inla.scale.model.internal(Q, constr = list(A = matrix(1, nrow = 1, ncol = n), e = 0))
    res <- gmrfs_scale_model_internal(Q, constr = list(A = matrix(1, nrow = 1, ncol = n), e = 0))
    Q <- res$Q
    f <- mean(res$var, na.rm = TRUE) - 1.0
    use.eigenvalues <- FALSE
  } else {
    n <- length(eigenvalues)
    eigenvalues <- pmax(0, sort(eigenvalues, decreasing = TRUE))
    gamma.invm1 <- c(1 / eigenvalues[1:(n - rankdef)], rep(0, rankdef)) - 1.0
    f <- mean(marginal.variances) - 1.0
    use.eigenvalues <- TRUE
  }

  if (use.eigenvalues) {
    ## this is fast for low dimension where we can compute the
    ## eigenvalues
    phi.s <- 1 / (1 + exp(-seq(-15, 12, length.out = 1000)))
    d <- numeric(length(phi.s))
    k <- 1
    for (phi in phi.s) {
      aa <- n * phi * f
      ## bb = sum(log(1-phi + phi*gamma.inv))
      bb <- sum(log1p(phi * gamma.invm1))
      ## bb = sum(log(1+phi*gamma.invm1))
      ## sometimes we *might* experience numerical problems, so we need
      ## to remove (if any) small phi's.
      d[k] <- (if (aa >= bb) sqrt(aa - bb) else NA)
      k <- k + 1
    }
  } else {
    ## alternative strategy for larger matrices
    phi.s <- 1 / (1 + exp(-c(seq(-15, 0, length.out = 40), 1:12)))
    d <- numeric(length(phi.s))
    log.q1.det <- gmrfs_sparse_det_bym(Q,
      adjust.for.con.comp = adjust.for.con.comp,
      constr = res$constr, rankdef = rankdef
    )
    d <- c(unlist(
      gmrfs_mclapply(
        phi.s,
        function(phi) {
          aa <- n * phi * f
          bb <- (n * log((1.0 - phi) / phi) +
            gmrfs_sparse_det_bym(Q + phi / (1 - phi) * Matrix::Diagonal(n),
              adjust.for.con.comp = adjust.for.con.comp,
              constr = res$constr, rankdef = rankdef
            ) -
            (log.q1.det - n * log(phi)))
          my.debug("aa=", aa, " bb=", bb, " sqrt(", aa - bb, ")")
          return(if (aa >= bb) sqrt(aa - bb) else NA)
        }
      )
    ))
  }
  names(d) <- NULL

  ## remove phi's that failed, if any
  remove <- is.na(d)
  my.debug(paste(sum(remove), "out of", length(remove), "removed"))
  d <- d[!remove]
  phi.s <- phi.s[!remove]
  remove <- 1:6
  d <- d[-remove]
  phi.s <- phi.s[-remove]
  phi.intern <- log(phi.s / (1 - phi.s))
  ff.d <- stats::splinefun(phi.intern, log(d))
  phi.intern <- seq(min(phi.intern) - 0, max(phi.intern), length.out = 10000)
  phi.s <- 1 / (1 + exp(-phi.intern))
  d <- exp(ff.d(phi.intern))

  ## ff.d:  return distance as a function of phi.intern
  ff.d.core <- stats::splinefun(phi.intern, log(d))
  ff.d <- function(phi.intern, deriv = 0) {
    if (deriv == 0) {
      return(exp(ff.d.core(phi.intern)))
    } else if (deriv == 1) {
      return(exp(ff.d.core(phi.intern)) * ff.d.core(phi.intern, deriv = 1))
    } else {
      stop("ERROR")
    }
  }
  ## f.d: return distance as a function of phi
  f.d <- function(phi.s) ff.d(log(phi.s / (1 - phi.s)))
  d <- f.d(phi.s)
  if (missing(lambda)) {
    ## Prob(phi < u) = alpha gives an analytical solution
    stopifnot(alpha > 0.0 && alpha < 1.0 && u > 0)
    lambda <- -log(1 - alpha) / f.d(u)
  }

  ## this is the log.prior wrt phi.intern
  log.jac <- log(abs(ff.d(phi.intern, deriv = 1)))
  log.prior <- log(lambda) - lambda * d + log.jac
  ssum <- sum(exp(log.prior) * (c(0, diff(phi.intern)) + c(diff(phi.intern), 0))) / 2.0
  my.debug("empirical integral of the prior = ", ssum, ". correct it.")
  log.prior <- log.prior - log(ssum)

  if (return.as.table) {
    ## return this prior as a "table:" converted to the intern scale for phi.
    my.debug("return prior as a table on phi.intern")
    theta.prior.table <- paste(c("table:", cbind(phi.intern, log.prior)),
      sep = "", collapse = " "
    )
    ## gmrfs: what the exact right tail needs beyond the table. For large
    ## logit(phi), d^2 grows by rankdef * log(1 + e^theta) (the null space of
    ## the Besag part), the other terms being O(e^-theta).
    attr(theta.prior.table, "tail") <- c(
      theta_max = phi.intern[length(phi.intern)], d_max = d[length(d)],
      rankdef = rankdef, lambda = lambda, log_ssum = log(ssum)
    )
    return(theta.prior.table)
  } else {
    ## return a function which evaluates the log-prior as a function of phi.
    my.debug("return prior as a function of phi")
    ## do the interpolation in the intern scale
    prior.theta.1 <- stats::splinefun(phi.intern, log.prior + phi.intern + 2 * log(1 + exp(-phi.intern)))
    prior.phi <- function(phi) prior.theta.1(log(phi / (1.0 - phi)))
    return(prior.phi)
  }
}
