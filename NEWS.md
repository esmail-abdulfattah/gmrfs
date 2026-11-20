# gmrfs 0.1.0

First CRAN release.

## Models

* `fbesag()`: Besag model with a separate precision per partition of the
  graph.
* `fbym2()`: BYM2 model whose spatial part is `fbesag()`.
* `cblocks()`, `pblocks()`, `sblocks()`: collapsed, partial and stacked
  block models, with shrinkage priors on the block scales.
* `omp_select()`: orthogonal matching pursuit screen for active blocks.

Each model can be used as `f(idx, model = .)` in INLA, or evaluated from R
with `gmrf_Q()`, `gmrf_mu()`, `gmrf_initial()`, `gmrf_log_prior()` and
`gmrf_log_norm_const()`. No model needs INLA to be built or evaluated:
`fbesag()` and `fbym2()` use the package's own copies of INLA's graph
reading, scaling and prior of the mixing parameter (with the marginal
variances by the Takahashi recursions in C). INLA is needed only to fit with
`INLA::inla()`.

## Use outside INLA (MCMC)

* `gmrf_log_density()`: `log p(x | theta)` for every model, with its linear
  constraints, computed as INLA does.
* `gmrf_sample()`: draws of `x | theta` on the constraints; with `z`, a
  non-centred parametrisation.
* Both take `Q_add` and `b` for `x | theta, y` with Gaussian data: exact
  Gibbs steps and the marginal likelihood with `x` integrated out.
* `demo("fbym2_mcmc")`: `fbym2` fitted by MCMC, compared with INLA.
* `gmrf_log_prior()` returns the log prior density of `theta`. For
  `fbesag()` and `fbym2()`, `inla = TRUE` gives the value passed to INLA
  instead, which adds a constant that puts INLA's marginal likelihood on the
  convention of `besag`/`bym2`.

## `fbym2()`: the prior of the mixing parameter

* The PC prior of `phi` is tabulated for `logit(phi)` up to 12, with the same
  values as INLA's `bym2`. Beyond, it now has the exact tail of the PC prior;
  it was extended linearly, which made it rise without bound, and INLA's
  `bym2` keeps it constant there. The marginal likelihood is unchanged when
  the posterior of `phi` stays inside the table.
* `gmrf_log_prior()` returns the prior of `phi` normalised on the whole line.
* `Q` stays finite for `phi` close to 1.

## `fbesag()` and `fbym2()`: marginal likelihood on large graphs

* INLA adds `diagonal = 1e-4` to the precision of a model with constraints.
  INLA computes the normalising constant of these models from that matrix,
  so the diagonal entered the marginal likelihood: 0.26 too high on a
  150 x 150 lattice. Both models now set `diagonal = 1e-8` in their handle
  (a `diagonal` given in `f()` still takes precedence), and their own jitter
  is `1e-8` instead of `1e-5`. On that lattice the marginal likelihood is now
  exact to 1e-3, where INLA's own `besag` is 0.26 too low with
  `diagonal = 1e-8` and 1.17 too low by default.
* `gmrf_sample()` applies the constraint correction twice, which removes the
  rounding error of a draw that is large along the null space.

## Changes to `fbesag()` since the GitHub version

* Graphs with several connected components and islands: one scaling
  constant and one sum-to-zero constraint per component, as INLA's `besag`;
  an island has precision equal to its partition's precision.
* The prior of the hyperparameters is a proper density, and with one
  partition the marginal likelihood is the same as INLA's `besag`. Use
  `nrep` with `f(..., replicate = )`.
* `P = 1` partition is allowed.
* Neighbours read from a graph file in any order are handled correctly.
* The C entry point is now `inla_cgeneric_gmrfs_fbesag`, so the package's
  code is used rather than the older copy built into INLA.
