# gmrfs

Adaptive Gaussian Markov random field models. Each model can be evaluated
directly from R (precision matrix, mean, hyperparameter prior, normalising
constant) for use in your own code, or used as a latent component in
[R-INLA](https://www.r-inla.org/) via the `cgeneric` interface. Bundles five
latent models:

| wrapper     | C entry point                 |
|-------------|-------------------------------|
| `fbesag()`  | `inla_cgeneric_gmrfs_fbesag`  |
| `fbym2()`   | `inla_cgeneric_gmrfs_fbym2`   |
| `cblocks()` | `inla_cgeneric_gmrfs_cblocks` |
| `pblocks()` | `inla_cgeneric_gmrfs_pblocks` |
| `sblocks()` | `inla_cgeneric_gmrfs_sblocks` |

## Use without INLA

```r
library(gmrfs)
cm <- cblocks(n = n, p = p, G_list = G_list)

theta <- gmrf_initial(cm)            # hyperparameters, internal scale
Q     <- gmrf_Q(cm, theta)           # sparse precision matrix
gmrf_mu(cm, theta)                   # mean
gmrf_log_prior(cm, theta)            # log prior of theta
gmrf_log_norm_const(cm, theta)       # log normalising constant of x | theta

x <- gmrf_sample(cm, theta, n = 10)  # draws of x | theta, constraints included
gmrf_log_density(cm, x, theta)       # log p(x | theta), constraints included
```

With these, the models can be used in your own MCMC sampler or optimiser.
`gmrf_sample()` and `gmrf_log_density()` also take `Q_add` and `b` for
`x | theta, y` with Gaussian data: exact Gibbs steps for `x`, and the
marginal likelihood `p(y | theta)` with `x` integrated out. See
`demo("fbym2_mcmc", package = "gmrfs")`, which fits `fbym2` by MCMC and
compares with INLA.

None of this needs INLA: every model, including `fbesag()` and `fbym2()`,
is built, evaluated and sampled without it. For those two the package has
its own copies of INLA's graph reading, scaling (`inla.scale.model`, with
the marginal variances by the same Takahashi recursions) and the prior of
the mixing parameter. INLA is needed only to fit with `INLA::inla()`.

## Install

```r
# INLA (if not already installed)
install.packages("INLA",
                 repos = c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                 dep = TRUE)

# devtools / remotes
install.packages("remotes")
remotes::install_github("esmail-abdulfattah/gmrfs")
```

The user's R will compile `src/*.c` at install time. A working `gcc` (with
OpenMP) plus BLAS/LAPACK is required, same toolchain expectations as any
package with C in `src/`.

## Quick start (collapsed)

```r
library(INLA); library(gmrfs)
# G_list: list of p covariance matrices, each n x n.
cm  <- cblocks(n = n, p = p, G_list = G_list, tau0 = 1.0)
fit <- inla(y ~ 1 + f(idx, model = cm),
            data = data.frame(y = y, idx = ind))
```

`tau0` is a horseshoe-style global scale on each block's standard deviation:
larger -> looser prior, smaller -> more aggressive shrinkage of inactive
blocks. A common heuristic is `tau0 = s / p` where `s` is the prior guess at
how many blocks are active.

See [`inst/examples/fit_collapsed.R`](inst/examples/fit_collapsed.R) for a
self-contained simulated example, and the `gmrfs_test/` companion repo for
the three formulations side-by-side.

## Try it (smoke tests)

Each model ships a tiny demo (`< 3 s`) that builds the cgeneric, runs INLA
end-to-end, and verifies `m$ok`:

```r
library(gmrfs)
demo(package = "gmrfs")            # list available demos
demo("fbesag",  package = "gmrfs") # 5-node partitioned Besag
demo("cblocks", package = "gmrfs") # collapsed blocks
demo("pblocks", package = "gmrfs") # partial blocks
demo("sblocks", package = "gmrfs") # stacked blocks
demo("fbym2_mcmc", package = "gmrfs") # fbym2 by MCMC vs INLA (~40 s)
```

A green `[OK] ... ran successfully.` line at the end of each confirms the
wrapper -> shlib -> C symbol -> INLA path is wired correctly.

## Models, input contracts

All three block models take a length-`p` list of per-block matrices at the
R layer; the wrappers handle stacking and any required factorisation.

- **`cblocks()`** (collapsed). Latent dim `n`. Covariance
  `Sigma(theta) = sum_j exp(-theta_j) G_j`.
  Inputs: `n`, `p`, `G_list` (length-`p` list of `n x n` covariances),
  `tau0`. `p` hyperparameters; `theta[j]` -> `G_list[[j]]` (no reorder).

- **`pblocks()`** (partial: one block separated, rest collapsed). Latent
  dim `2 * n`. Block `separate` (default 1) lives at indices `1..n` with
  precision `exp(theta[1]) * R_separate` (where `R_separate =
  G_list[[separate]]^{-1}`); the remaining `p - 1` blocks live at
  `n+1..2n` with covariance `sum_{j != separate} exp(-theta[j]) G_j`.
  Inputs: `n`, `p`, `G_list`, `separate`, `tau0`. `p` hyperparameters; the
  mapping back to `G_list` is recorded in `attr(cm, "block_order")`.

- **`sblocks()`** (stacked). Latent dim `p * n`. Precision is
  block-diagonal, `diag(exp(theta_1) R_1, ..., exp(theta_p) R_p)`.
  Inputs: `n`, `p`, `R_list` (length-`p` list of `n x n` precisions,
  typically `R_j = G_j^{-1}`), `tau0`. `p` hyperparameters.

- **`fbesag()`** (partitioned Besag with `P` partitions). Inputs:
  `graph`, `id` (region -> partition map), `sd_gamma`, PC-prior
  `param = list(p1, p2)`. `P` hyperparameters.

## License

GPL (>= 2). Includes code from R-INLA (MIT licence); see
[`inst/COPYRIGHTS`](inst/COPYRIGHTS).
