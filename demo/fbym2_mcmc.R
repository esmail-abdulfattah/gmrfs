## fbym2 fitted by MCMC, without INLA's inference: only gmrf_log_density(),
## gmrf_sample() and gmrf_log_prior(). INLA is used to build the model and,
## at the end, to compare with its posterior.
##
##   y_ij ~ N(v_i, 1/tau_y), j = 1..4,  x = (v, u) ~ fbym2(theta),
##   theta = (log tau, logit phi, d_1), tau_y ~ Gamma(1, 5e-5) (INLA's default)
##
## The field is integrated out exactly. With Q_add = tau_y B'B and
## b = tau_y B'y (B maps x to the mean of each observation), at any x on the
## constraints
##   log p(y | theta) = log p(x | theta) + log p(y | x) - log p(x | y, theta)
## and x | y, theta is drawn exactly. So the sampler is
##   (log tau_y, theta)  random-walk Metropolis on log p(y | theta) + log priors
##   x | y, theta        exact draws, for a thinned subset of the chain

suppressPackageStartupMessages({
  library(gmrfs)
  library(Matrix)
  library(INLA)
})
options(mc.cores = 1)
set.seed(1)

## ---- 10 x 10 lattice, two partitions (left half / right half) -------------
nr <- 10; n <- nr * nr
W <- matrix(0, n, n)
for (i in 1:nr) for (j in 1:nr) {
  k <- (i - 1) * nr + j
  if (j < nr) W[k, k + 1] <- W[k + 1, k] <- 1
  if (i < nr) W[k, k + nr] <- W[k + nr, k] <- 1
}
id <- ifelse(rep(1:nr, times = nr) <= nr / 2, 1, 2)
m <- fbym2(W, id, sd_gamma = 1)
N <- 2 * n
n_obs <- 4                                   # observations per area: the noise
obs <- rep(1:n, each = n_obs)                # is then told apart from bym2's iid part
B <- sparseMatrix(i = seq_along(obs), j = obs, x = 1, dims = c(length(obs), N))

## ---- data, simulated from the model itself ---------------------------------
theta_true <- c(log_tau_y = log(4), log_tau = 0, logit_phi = qlogis(0.8), d_1 = -1)
x_true <- gmrf_sample(m, theta_true[-1])[, 1]
y <- rnorm(length(obs), x_true[obs], 1 / sqrt(exp(theta_true[1])))

## ---- log posterior of (log tau_y, theta), x integrated out -----------------
log_post <- function(th) {
  tau_y <- exp(th[1]); theta <- th[-1]
  Q_add <- tau_y * crossprod(B); b <- tau_y * as.numeric(crossprod(B, y))
  ## the posterior mean of x (z = 0), a point on the constraints
  x <- gmrf_sample(m, theta, z = numeric(N), Q_add = Q_add, b = b)[, 1]
  gmrf_log_density(m, x, theta) +
    sum(stats::dnorm(y, x[obs], 1 / sqrt(tau_y), log = TRUE)) -
    gmrf_log_density(m, x, theta, Q_add = Q_add, b = b) +
    stats::dgamma(tau_y, 1, 5e-5, log = TRUE) + th[1] +        # prior of log tau_y
    gmrf_log_prior(m, theta)
}

## ---- random-walk Metropolis, step sizes adapted during burn-in -------------
n_iter <- 8000; burn <- 3000
th <- c(0, gmrf_initial(m)); th[2:3] <- 0
lp <- log_post(th)
step <- rep(0.3, 4); win <- 0
keep <- matrix(NA, n_iter, 4, dimnames = list(NULL, names(theta_true)))
t0 <- Sys.time()
for (it in 1:n_iter) {
  prop <- th + stats::rnorm(4, 0, step)
  lp_prop <- log_post(prop)
  if (log(stats::runif(1)) < lp_prop - lp) { th <- prop; lp <- lp_prop; win <- win + 1 }
  keep[it, ] <- th
  if (it <= burn && it %% 100 == 0) {                # aim at acceptance 0.25
    step <- step * exp(win / 100 - 0.25); win <- 0
  }
}
post <- keep[-(1:burn), ]
cat(sprintf("MCMC: %d iterations in %.0f s, acceptance after burn-in %.2f\n", n_iter,
            as.numeric(Sys.time() - t0, units = "secs"),
            mean(diff(post[, 1]) != 0)))

## x | y, theta for every 10th draw
idx <- seq(1, nrow(post), by = 10)
v <- sapply(idx, function(k) {
  tau_y <- exp(post[k, 1])
  gmrf_sample(m, post[k, -1], Q_add = tau_y * crossprod(B), b = tau_y * as.numeric(crossprod(B, y)))[1:n, 1]
})

## ---- compare with INLA -----------------------------------------------------
fit <- inla(y ~ -1 + f(idx, model = m), data = list(y = y, idx = obs))
h <- fit$internal.summary.hyperpar
tab <- data.frame(truth = theta_true,
                  MCMC_mean = colMeans(post), INLA_mean = h$mean,
                  MCMC_sd = apply(post, 2, stats::sd), INLA_sd = h$sd)
print(round(tab, 3))
cat(sprintf("latent v: max |MCMC mean - INLA mean| = %.3f (posterior sd ~ %.2f)\n",
            max(abs(rowMeans(v) - fit$summary.random$idx$mean[1:n])),
            mean(fit$summary.random$idx$sd[1:n])))
