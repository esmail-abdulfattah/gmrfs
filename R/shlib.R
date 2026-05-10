## Resolve which shared library to hand to inla.cgeneric.define.
##
## Priority (highest first):
##   1. explicit `shlib` argument (always wins).
##   2. INLA's precompiled external library: INLA::inla.external.lib(name).
##      Only consulted when `prefer_inla = TRUE` (the default); when INLA
##      ships a copy, that is the canonical one.
##   3. local-dev directory: getOption("gmrfs.local_dir") or env var
##      GMRFS_LOCAL_DIR.  If set, look for <name><dynlib.ext> there
##      (e.g. /path/to/gmrfs/src/cblocks.so).  Used when `prefer_inla`
##      is FALSE, or when INLA does not (yet) ship the model.
##   4. the package's own shared library compiled at install time.

gmrfs_resolve_shlib <- function(name, shlib = NULL, prefer_inla = TRUE) {

  ## (1) explicit per-call override
  if (!is.null(shlib)) {
    if (!is.character(shlib) || length(shlib) != 1L || !nzchar(shlib) ||
        !file.exists(shlib)) {
      stop(sprintf("explicit shlib not found: %s", shlib), call. = FALSE)
    }
    return(normalizePath(shlib))
  }

  ## (2) INLA's precompiled external library (only if prefer_inla)
  if (isTRUE(prefer_inla)) {
    p <- tryCatch(INLA::inla.external.lib(name), error = function(e) NULL)
    if (!is.null(p) && nzchar(p) && file.exists(p)) return(p)
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

gmrfs_own_shlib <- function() {
  dll <- getLoadedDLLs()[["gmrfs"]]
  if (is.null(dll)) {
    stop("gmrfs shared library is not loaded; reinstall the package",
         call. = FALSE)
  }
  dll[["path"]]
}
