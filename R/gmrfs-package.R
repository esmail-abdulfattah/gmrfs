#' gmrfs: Adaptive GMRF models for INLA via cgeneric
#'
#' Provides four cgeneric latent models — `fbesag`,
#' `cblocks`, `pblocks`, `sblocks` — for use as
#' `f(idx, model = .)` in [INLA::inla()] formulas. When INLA ships a
#' precompiled library for a model (queried via
#' [INLA::inla.external.lib()]), the wrappers use that copy; otherwise
#' they use the shared library compiled at install time from
#' `src/*.c`.
#'
#' @name gmrfs-package
#' @aliases gmrfs
#' @keywords internal
"_PACKAGE"
