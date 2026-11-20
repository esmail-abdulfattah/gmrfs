## Graphs and the scaling of intrinsic models, without INLA.
##
## Copied from R-INLA (rinla/R/read.graph.R, graph.matrix.R, sm.R and
## scale.model.R), renamed: inla.read.graph(), inla.add.graph.cc(),
## inla.graph2matrix(), inla.as.sparse() and inla.scale.model.internal().
## R-INLA is Copyright
## Haavard Rue, MIT licence. Changed: the marginal variances, which INLA gets
## from its binary (inla.qinv()), are computed by gmrfs_qinv_diag() below, with
## the same algorithm (Takahashi recursions, src/qinv.c); spdep "nb" objects
## are read without spdep.

## inla.as.sparse(): a general triplet matrix (dgTMatrix)
gmrfs_as_sparse <- function(A, unique = TRUE) {
  if (!inherits(A, "Matrix")) A <- methods::as(A, "Matrix")
  if (unique) {
    A <- methods::as(methods::as(methods::as(methods::as(A, "dMatrix"), "generalMatrix"),
                                 "CsparseMatrix"), "TsparseMatrix")
  } else if (!methods::is(A, "dgTMatrix")) {
    A <- methods::as(methods::as(methods::as(A, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  }
  A
}

## inla.graph.binary.file.magic()
.gmrfs_graph_magic <- -1L

## inla.read.graph(): a graph file (ascii or binary), the lines of a graph
## file, a (sparse) matrix whose nonzero off-diagonal entries are the edges,
## an "inla.graph", or an spdep "nb" object. Returns an "inla.graph": n,
## nnbs, nbs, and the connected components cc.
gmrfs_read_graph <- function(..., size.only = FALSE) {
  read_ascii <- function(filename, offset = 0L, size.only = FALSE) {
    stopifnot(file.exists(filename))
    s <- readLines(filename)
    if (length(s) == 0L) return(NULL)
    s <- sapply(s, function(x) gsub("#.*$", "", x))
    s <- as.integer(unlist(sapply(s, function(x) strsplit(x, "[ \t]+"))))
    s <- s[!is.na(s)]
    n <- s[1L]
    if (size.only) return(n)
    g <- list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n))
    k <- 2L
    for (i in 1L:n) {
      ## a graph file may number the nodes from 0
      if (s[k] + offset == 0L) return(read_ascii(filename, offset = 1L))
      stopifnot(s[k] + offset >= 1L && s[k] + offset <= n)
      idx <- s[k] + offset
      k <- k + 1L
      g$nnbs[idx] <- s[k]
      k <- k + 1L
      if (g$nnbs[idx] > 0L) {
        g$nbs[[idx]] <- s[k:(k + g$nnbs[idx] - 1L)] + offset
        k <- k + g$nnbs[idx]
      }
    }
    stopifnot(k - 1L == length(s))
    class(g) <- "inla.graph"
    gmrfs_add_graph_cc(g)
  }

  read_binary <- function(filename, offset = 0L, size.only = FALSE) {
    stopifnot(file.exists(filename))
    fp <- gzfile(filename, "rb")
    s <- as.integer(readBin(fp, integer(), n = 1L))
    close(fp)
    if (length(s) == 0L || s[1L] != .gmrfs_graph_magic) return(NULL)
    n.try <- 2^12
    repeat {
      fp <- gzfile(filename, "rb")
      s <- as.integer(readBin(fp, integer(), n = n.try))
      close(fp)
      if (length(s) < n.try) break else n.try <- n.try * 4L
    }
    s <- s[-1L]
    n <- s[1L]
    if (size.only) return(n)
    g <- list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n))
    k <- 2L
    for (i in 1L:n) {
      if (s[k] + offset == 0L) return(read_binary(filename, offset = 1L))
      stopifnot(s[k] + offset >= 1L && s[k] + offset <= n)
      idx <- s[k] + offset
      k <- k + 1L
      g$nnbs[idx] <- s[k]
      k <- k + 1L
      if (g$nnbs[idx] > 0L) {
        g$nbs[[idx]] <- s[k:(k + g$nnbs[idx] - 1L)] + offset
        k <- k + g$nnbs[idx]
      }
    }
    stopifnot(k - 1L == length(s))
    class(g) <- "inla.graph"
    gmrfs_add_graph_cc(g)
  }

  matrix2graph <- function(Q, size.only = FALSE) {
    n <- dim(Q)
    if (size.only) return(n[1L])
    if (n[1] != n[2]) stop(paste("Matrix must be a square matrix, dim(Q) =", dim(Q)))
    n <- n[1]
    g <- list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n), graph.file = NA)
    Q <- gmrfs_as_sparse(Q)
    diag(Q) <- 1
    Q <- gmrfs_as_sparse(Q)
    ord <- order(Q@i)
    Q@i <- Q@i[ord]; Q@j <- Q@j[ord]; Q@x <- Q@x[ord]
    hash.len <- table(Q@i)
    hash.idx <- c(1L, 1L + cumsum(hash.len))
    stopifnot(length(hash.len) == ncol(Q))
    for (i in 1L:n) {
      if (hash.len[i] > 1L) {
        idx <- hash.idx[i]:(hash.idx[i] + hash.len[i] - 1L)
        j <- Q@j[idx] + 1L
        x <- Q@x[idx]
        j <- j[(x != 0) & (j != i)]
      } else {
        j <- numeric()
      }
      ## (INLA assigns NULL here, which drops the list element and relies on
      ## padding the list afterwards; keep an empty element instead)
      g$nbs[i] <- list(j)
      g$nnbs[i] <- length(j)
    }
    class(g) <- "inla.graph"
    gmrfs_add_graph_cc(g)
  }

  args <- list(...)
  graph <- args[[1L]]
  if (is.character(graph) || length(args) > 1L ||
      (is.numeric(graph) && !(is.matrix(graph) || methods::is(graph, "Matrix")))) {
    graph <- paste(as.character(graph))
    if (length(graph) == 1L && file.exists(graph)) {
      g <- read_binary(graph, size.only = size.only)
      if (is.null(g)) g <- read_ascii(graph, size.only = size.only)
      return(g)
    }
    tfile <- tempfile()
    on.exit(unlink(tfile))
    cat(unlist(args), sep = "\n", file = tfile, append = FALSE)
    return(gmrfs_read_graph(tfile, size.only = size.only))
  }
  if (inherits(graph, "inla.graph")) {
    if (size.only) return(graph$n)
    if (is.null(graph$cc)) graph <- gmrfs_add_graph_cc(graph)
    return(graph)
  }
  if (inherits(graph, "nb")) {
    ## an spdep "nb": a list of neighbour indices, 0 for none
    n <- length(graph)
    if (size.only) return(n)
    i <- rep(seq_len(n), lengths(graph))
    j <- unlist(graph)
    keep <- j > 0
    W <- Matrix::sparseMatrix(i = i[keep], j = j[keep], x = 1, dims = c(n, n))
    return(matrix2graph(W + Matrix::t(W)))
  }
  matrix2graph(graph, size.only = size.only)
}

## inla.add.graph.cc(): the connected components of a graph
gmrfs_add_graph_cc <- function(graph) {
  cc <- list(id = NA, n = NA, nodes = NA)
  n <- graph$n
  s <- integer(n)
  do.visit <- function(idxs, k) {
    repeat {
      idxs.visit <- c()
      if (length(idxs) > 0) {
        for (idx in idxs) {
          if (s[idx] == 0L) {
            s[idx] <<- k
            visit.next <- graph$nbs[[idx]]
            idxs.visit <- c(idxs.visit, visit.next[which(s[visit.next] == 0L)])
          }
        }
        idxs <- unique(sort(idxs.visit))
      } else {
        break
      }
    }
    invisible()
  }
  k <- 1L
  for (i in 1L:n) {
    if (s[i] == 0L) {
      do.visit(i, k)
      k <- k + 1L
    }
  }
  cc$id <- s
  cc$n <- max(s)
  cc$nodes <- lapply(1L:cc$n, function(cc.id, cs) sort(which(cc.id == cs)), cs = s)
  for (ii in which(lengths(cc$nodes) == 1)) s[cc$nodes[[ii]]] <- NA
  cc$mean <- factor(s, exclude = NA)
  graph$cc <- cc
  graph
}

## inla.graph2matrix(): the pattern matrix of a graph, 1 on the diagonal
gmrfs_graph2matrix <- function(graph) {
  g <- gmrfs_read_graph(graph)
  i <- rep(1:g$n, 1L + g$nnbs)
  j <- unlist(lapply(1:g$n, function(i) c(i, g$nbs[[i]])))
  stopifnot(length(i) == length(j))
  gmrfs_as_sparse(Matrix::sparseMatrix(i = i, j = j, x = rep(1, length(i)), dims = c(g$n, g$n)))
}

## diag(Q^-1), or with constraints A x = e, diag of the covariance of
## x | A x = e: Sigma - Sigma A' (A Sigma A')^-1 A Sigma. As inla.qinv(), from
## the Cholesky factor by the Takahashi recursions (src/qinv.c).
gmrfs_qinv_diag <- function(Q, constr = NULL) {
  Q <- Matrix::forceSymmetric(methods::as(methods::as(Q, "dMatrix"), "CsparseMatrix"), uplo = "U")
  ch <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE, super = FALSE)
  d <- numeric(nrow(Q))
  d[ch@perm + 1L] <- .Call(gmrfs_takahashi_diag, ch@p, ch@i, ch@x, ch@nz)
  if (!is.null(constr) && NROW(constr$A) > 0) {
    A <- matrix(constr$A, ncol = nrow(Q))
    W <- as.matrix(Matrix::solve(ch, t(A)))
    d <- d - rowSums((W %*% solve(A %*% W)) * W)
  }
  d
}

## inla.scale.model.internal(): scale each connected component of Q so that
## the geometric mean of its marginal variances (under the constraints that
## touch it) is 1; a component of size 1 gets Q = 1. Returns the scaled Q and
## the scaled marginal variances.
gmrfs_scale_model_internal <- function(Q, constr = NULL, eps = sqrt(.Machine$double.eps)) {
  marg.var <- rep(0, nrow(Q))
  Q <- gmrfs_as_sparse(Q)
  g <- gmrfs_read_graph(Q)
  for (k in seq_len(g$cc$n)) {
    i <- g$cc$nodes[[k]]
    n <- length(i)
    QQ <- Q[i, i, drop = FALSE]
    if (n == 1) {
      QQ[1, 1] <- 1
      marg.var[i] <- 1
    } else {
      cconstr <- constr
      if (!is.null(constr)) {
        cconstr$A <- constr$A[, i, drop = FALSE]
        eeps <- eps
        idx.zero <- which(rowSums(abs(cconstr$A)) == 0)
        if (length(idx.zero) > 0) {
          cconstr$A <- cconstr$A[-idx.zero, , drop = FALSE]
          cconstr$e <- cconstr$e[-idx.zero]
        }
      } else {
        eeps <- 0
      }
      res <- gmrfs_qinv_diag(QQ + Matrix::Diagonal(n) * max(Matrix::diag(QQ)) * eeps, constr = cconstr)
      fac <- exp(mean(log(res)))
      QQ <- fac * QQ
      marg.var[i] <- res / fac
    }
    Q[i, i] <- QQ
  }
  list(Q = Q, var = marg.var)
}

## inla.scale.model()
gmrfs_scale_model <- function(Q, constr = NULL, eps = sqrt(.Machine$double.eps)) {
  gmrfs_scale_model_internal(Q = Q, constr = constr, eps = eps)$Q
}
