suppressPackageStartupMessages({
  library(gmrfs); library(INLA)
})

cm <- fbesag(graph    = matrix(rep(1, 5), 5, 5),
             id       = c(1, 1, 1, 2, 2),
             sd_gamma = 0.15,
             param    = list(p1 = 1, p2 = 1e-5))

set.seed(1)
m <- INLA::inla(y ~ 1 + f(idx, model = cm),
                data   = list(y = rnorm(5, 0, 0.1), idx = 1:5),
                family = "gaussian")

print(summary(m)$fixed)
stopifnot(isTRUE(m$ok))
cat("\n\n[OK] fbesag minimal demo ran successfully.\n")
