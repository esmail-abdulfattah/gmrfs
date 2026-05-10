# gmrfs

Adaptive GMRF models for [R-INLA](https://www.r-inla.org/) via the `cgeneric`
interface. Bundles four latent models:

| wrapper     | C entry point            | shipped in INLA? |
|-------------|--------------------------|------------------|
| `fbesag()`  | `inla_cgeneric_fbesag`   | yes (today)      |
| `cblocks()` | `inla_cgeneric_cblocks`  | not yet — built locally |
| `pblocks()` | `inla_cgeneric_pblocks`  | not yet — built locally |
| `sblocks()` | `inla_cgeneric_sblocks`  | not yet — built locally |

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
OpenMP) plus BLAS/LAPACK is required — same toolchain expectations as any
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
```

A green `[OK] ... ran successfully.` line at the end of each confirms the
wrapper -> shlib -> C symbol -> INLA path is wired correctly.

## Models — input contracts

All three block models take a length-`p` list of per-block matrices at the
R layer; the wrappers handle stacking and any required factorisation.

- **`cblocks()`** (collapsed). Latent dim `n`. Covariance
  `Sigma(theta) = sum_j exp(-theta_j) G_j`.
  Inputs: `n`, `p`, `G_list` (length-`p` list of `n x n` covariances),
  `tau0`. `p` hyperparameters; `theta[j]` -> `G_list[[j]]` (no reorder).

- **`pblocks()`** (partial — one block separated, rest collapsed). Latent
  dim `2 * n`. Block `separate` (default 1) lives at indices `1..n` with
  precision `exp(theta[1]) * R_separate` (where `R_separate =
  G_list[[separate]]^{-1}`); the remaining `p - 1` blocks live at
  `n+1..2n` with covariance `sum_{j != separate} exp(-theta[j]) G_j`.
  Inputs: `n`, `p`, `G_list`, `separate`, `tau0`. `p` hyperparameters; the
  mapping back to `G_list` is recorded in `attr(cm, "block_order")`.

- **`sblocks()`** (stacked). Latent dim `p * n`. Precision is
  block-diagonal, `diag(exp(theta_1) R_1, ..., exp(theta_p) R_p)`.
  Inputs: `n`, `p`, `R_list` (length-`p` list of `n x n` precisions —
  typically `R_j = G_j^{-1}`), `tau0`. `p` hyperparameters.

- **`fbesag()`** (partitioned Besag with `P` partitions). Inputs:
  `graph`, `id` (region -> partition map), `sd_gamma`, PC-prior
  `param = list(p1, p2)`. `P` hyperparameters.

## Identifiability — check before fitting

Variance components from a sum `s = sum_j u_j` with `u_j ~ N(0,
sigma_j^2 G_j^{-1})` are only jointly identifiable if the `G_j` are
*shaped differently*. If two `G_j` have similar eigenstructure, the
collapsed model will smear variance across them.

The Frobenius cosine matrix between blocks is the right pre-flight check:

```r
C_jk <- function(A, B) sum(A * B) / sqrt(sum(A * A) * sum(B * B))
C    <- outer(seq_along(G_list), seq_along(G_list),
              Vectorize(function(j, k) C_jk(G_list[[j]], G_list[[k]])))
e    <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
```

`C[j,j] = 1`; off-diagonals near 1 mean blocks `j` and `k` look the same.
Smallest eigenvalue of `C` (or its condition number) is the single number
to watch:

| smallest eigenvalue | verdict      |
|---------------------|--------------|
| > 0.1               | OK           |
| 0.01 -- 0.1         | MARGINAL     |
| < 0.01              | CONFOUNDED — drop or restructure a block |

The `gmrfs_test/` companion repo bundles a `check_block_identifiability()`
helper that prints `C`, its eigenvalues, and per-block effective rank.

## License

GPL (>= 2). Includes `src/cgeneric.h` from R-INLA.
