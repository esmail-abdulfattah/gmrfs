#' gmrfs: Adaptive Gaussian Markov random field models
#'
#' Provides five latent Gaussian Markov random field models: `fbesag`,
#' `fbym2`, `cblocks`, `pblocks` and `sblocks`.
#'
#' Each model can be used in two ways:
#'
#' * directly from R, through [gmrf_Q()], [gmrf_mu()], [gmrf_log_prior()]
#'   and the other functions in [gmrf_eval], for use in your own sampler
#'   or optimiser. This does not need INLA.
#' * as `f(idx, model = .)` in an `INLA::inla()` formula, through INLA's
#'   `cgeneric` interface. When INLA ships a precompiled library for a
#'   model, the wrappers use that copy; otherwise they use the shared
#'   library compiled with this package.
#'
#' None of the models needs INLA to be built or evaluated. For `fbesag()` and
#' `fbym2()` the package has its own copies of INLA's graph reading, scaling
#' and the prior of the mixing parameter; INLA is needed only to fit with
#' `INLA::inla()`.
#'
#' @name gmrfs-package
#' @aliases gmrfs
#' @keywords internal
"_PACKAGE"
