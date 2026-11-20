## Resolve which shared library to hand to inla.cgeneric.define.
##
## Priority (highest first):
##   1. explicit `shlib` argument (always wins).
##   2. INLA's precompiled external library: INLA::inla.external.lib(name).
##      Only consulted when `prefer_inla = TRUE` (the default), and only
##      used if it contains the function `symbol` that the model calls (a
##      copy under another name, or none, falls through to 3 and 4).
##   3. local-dev directory: getOption("gmrfs.local_dir") or env var
##      GMRFS_LOCAL_DIR.  If set, look for <name><dynlib.ext> there
##      (e.g. /path/to/gmrfs/src/cblocks.so).  Used when `prefer_inla`
##      is FALSE, or when INLA does not (yet) ship the model.
##   4. the package's own shared library compiled at install time.

gmrfs_resolve_shlib <- function(name, shlib = NULL, prefer_inla = TRUE, symbol = NULL) {

  ## (1) explicit per-call override
  if (!is.null(shlib)) {
    if (!is.character(shlib) || length(shlib) != 1L || !nzchar(shlib) ||
        !file.exists(shlib)) {
      stop(sprintf("explicit shlib not found: %s", shlib), call. = FALSE)
    }
    return(normalizePath(shlib))
  }

  ## (2) INLA's precompiled external library (only if prefer_inla)
  if (FALSE) {
    ## old: inla.external.lib() takes its argument's name, not its value
    ## (as library()), so this always looked for .../external/name/libname.so
    if (isTRUE(prefer_inla)) {
      p <- tryCatch(INLA::inla.external.lib(name), error = function(e) NULL)
      if (!is.null(p) && nzchar(p) && file.exists(p)) return(p)
    }
  }
  if (isTRUE(prefer_inla)) {
    p <- .gmrfs_inla_lib(name)
    if (!is.null(p) && (is.null(symbol) || .gmrfs_has_symbol(p, symbol))) return(p)
  }

  ## (3) local-dev directory
  local_dir <- getOption("gmrfs.local_dir",
                         default = Sys.getenv("GMRFS_LOCAL_DIR", unset = ""))
  if (nzchar(local_dir)) {
    candidate <- file.path(local_dir, paste0(name, .Platform$dynlib.ext))
    if (file.exists(candidate)) return(normalizePath(candidate))
  }

  ## (4) package's own .so/.dll built at install
  gmrfs_own_shlib()
}

## Path of INLA's precompiled library for `name`, or NULL if INLA is not
## installed or does not ship one. do.call() passes the value of `name`.
.gmrfs_inla_lib <- function(name) {
  if (!requireNamespace("INLA", quietly = TRUE)) return(NULL)
  p <- tryCatch(do.call(INLA::inla.external.lib, list(name)), error = function(e) NULL)
  if (is.null(p) || !nzchar(p) || !file.exists(p)) NULL else p
}

## Does the shared library at `path` export the function `symbol`? Loaded
## with local = TRUE and unloaded again, so nothing stays in the session.
.gmrfs_has_symbol <- function(path, symbol) {
  dll <- tryCatch(dyn.load(path, local = TRUE, now = FALSE), error = function(e) NULL)
  if (is.null(dll)) return(FALSE)
  on.exit(try(dyn.unload(path), silent = TRUE))
  is.loaded(symbol, PACKAGE = dll[["name"]])
}

## Map a prior name to the integer code understood by the cgeneric C side
## (priors.h). Integers are passed reliably to cgeneric models across INLA
## versions, unlike character arguments.
gmrfs_prior_code <- function(prior) {
  codes <- c(horseshoe = 1L, hs = 1L, pc.prec = 2L, pcprec = 2L,
             loggamma = 3L, gamma = 4L, gaussian = 5L, normal = 5L, flat = 6L)
  key <- tolower(prior)
  if (!key %in% names(codes))
    stop(sprintf("unknown prior '%s'; one of: %s", prior,
                 paste(unique(names(codes)), collapse = ", ")), call. = FALSE)
  unname(codes[key])
}

## Resolve prior parameters: if `param` is NULL, supply the default for the
## named prior. For the horseshoe (the block models' default) the single
## parameter is the global scale `tau0`; other priors require explicit `param`.
gmrfs_prior_param <- function(prior, param, tau0 = 1.0) {
  if (!is.null(param)) return(as.numeric(param))
  pl <- tolower(prior)
  if (pl %in% c("horseshoe", "hs", "flat")) return(as.numeric(tau0))
  if (pl == "loggamma") return(c(1, 0.00005))
  stop(sprintf("prior '%s' requires `param` (no default available)", prior),
       call. = FALSE)
}

gmrfs_own_shlib <- function() {
  dll <- getLoadedDLLs()[["gmrfs"]]
  if (is.null(dll)) {
    stop("gmrfs shared library is not loaded; reinstall the package",
         call. = FALSE)
  }
  dll[["path"]]
}
