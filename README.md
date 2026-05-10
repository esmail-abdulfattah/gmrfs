# gmrfs

Adaptive GMRF models for [R-INLA](https://www.r-inla.org/) via the `cgeneric`
interface. Bundles four latent models:

| wrapper | C entry point | shipped in INLA? |
|---|---|---|
| `fbesag()`    | `inla_cgeneric_fbesag`       | yes (today) |
| `cblocks()` | `inla_cgeneric_collapsed`    | not yet — built locally |
| `pblocks()`   | `inla_cgeneric_partial`      | not yet — built locally |
| `sblocks()`   | `inla_cgeneric_stacked`      | not yet — built locally |

For each call, the wrapper first asks `INLA::inla.external.lib(name)`; if INLA
ships a precompiled library for that model the wrapper uses it. Otherwise it
falls back to the shared library compiled at install time from `src/*.c`.
This means the package keeps working today and stays forward-compatible: as
the remaining three models get upstreamed into INLA, `prefer_inla = TRUE`
will pick up the shipped copy automatically without any wrapper change.

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
# G_list: list of n x n covariance matrices, one per block
cm <- cblocks(n = n, p = p, G = do.call(rbind, G_list))
fit <- inla(y ~ 1 + f(idx, model = cm),
            data = data.frame(y = y, idx = seq_len(n)))
```

See [`inst/examples/fit_collapsed.R`](inst/examples/fit_collapsed.R) for a
self-contained simulated example.

## Models — input contracts

- **collapsed**: `n` (block size, latent dim), `p` (#blocks), `G` is
  `(p*n) x n` row-stacked dense covariances. `p` hyperparameters (one
  log-precision per block).
- **partial**: latent dim `2*n`. Block 1 lives at indices `1..n` with
  precision `exp(theta_1) R1`; the rest live at `n+1..2n` with precision
  `(sum_{j>=2} exp(-theta_j) G_j)^{-1}`. Inputs: `n`, `p`, `G` `(p*n) x n`,
  sparse `R1` `n x n`. `p` hyperparameters.
- **stacked**: latent dim `p*n`. Precision is block-diagonal,
  `diag(exp(theta_j) R)` per block. Inputs: `n`, `p`, sparse `R` `(p*n) x
  (p*n)`. `p` hyperparameters.
- **fbesag**: partitioned Besag with `P` partitions.
  Inputs: `graph`, `id` (region→partition map), `sd_gamma`, PC-prior
  `param=list(p1, p2)`. `P` hyperparameters.

## License

GPL (>= 2). Includes `src/cgeneric.h` from R-INLA.
