#ifndef GMRFS_PRIORS_H
#define GMRFS_PRIORS_H

/*
 * priors.h -- a small catalog of hyperparameter prior log-densities for the
 * gmrfs cgeneric models, mirroring R-INLA's implementations:
 *   r-inla-devel/inlaprog/src/inla-priors.c   (loggamma, gamma, gaussian, ...)
 *   r-inla-devel/inlaprog/src/pc-priors.c     (pc.prec, ...)
 *
 * Convention: densities are on the INTERNAL scale of the hyperparameter.  For a
 * precision hyperparameter the internal representation is theta = log(prec),
 * which is what INLA passes; the PC and loggamma priors below assume that.
 *
 * Use gmrfs_log_prior(name, theta, params, nparams) to select a prior by name.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- individual priors (log density) ----------------------------------- */

/* Numerically stable e^x * E_1(x), x > 0 (used by the horseshoe prior). */
static inline double gmrfs_expE1(double x)
{
	if (!(x > 0.0)) return NAN;
	if (x < 1.0) {
		const double EG = 0.5772156649015329;
		double e1 = -EG - log(x);
		double term = 1.0;
		for (int k = 1; k < 100; k++) {
			term *= -x / (double) k;
			double add = -term / (double) k;
			e1 += add;
			if (fabs(add) < 1e-16 * (1.0 + fabs(e1))) break;
		}
		return exp(x) * e1;
	} else {
		const double TINY = 1e-30;
		double b = x + 1.0, d = 1.0 / b, c = 1.0 / TINY, h = d;
		for (int i = 1; i < 200; i++) {
			double a = -(double) (i * i);
			b += 2.0;
			double denom = b + a * d;
			if (fabs(denom) < TINY) denom = TINY;
			d = 1.0 / denom;
			double cnew = b + a / c;
			if (fabs(cnew) < TINY) cnew = TINY;
			c = cnew;
			double del = c * d;
			h *= del;
			if (fabs(del - 1.0) < 1e-15) break;
		}
		return h;
	}
}

/* Horseshoe-style half-Cauchy prior on a log-precision (the prior currently
 * used by cblocks/pblocks/sblocks): param = (tau0), the global scale. */
static inline double gmrfs_prior_horseshoe(double theta, double tau0)
{
	static const double LOG_PI_SQRT_2PI = 1.8378770664093453;
	double tau02 = tau0 * tau0;
	double x = (0.5 / tau02) * exp(-theta);
	double v = gmrfs_expE1(x);
	return -0.5 * theta + log(v) - LOG_PI_SQRT_2PI - 2.0 * log(tau0);
}

/* PC prior on a precision (INLA "pc.prec"): param = (u, alpha) with
 * P(sigma > u) = alpha,  sigma = 1/sqrt(prec) = exp(-theta/2).
 * lambda = -log(alpha)/u ;  log pi = log(lambda/2) - lambda*exp(-theta/2) - theta/2. */
static inline double gmrfs_prior_pc_prec(double theta, double u, double alpha)
{
	double lam = -log(alpha) / u;
	double t2 = 0.5 * theta;
	return log(0.5 * lam) - lam * exp(-t2) - t2;
}

/* loggamma on log-precision (INLA "loggamma"): prec = exp(theta) ~
 * Gamma(shape = a, rate = b).  log pi = a*log(b) - lgamma(a) + a*theta - b*exp(theta). */
static inline double gmrfs_prior_loggamma(double theta, double a, double b)
{
	return a * log(b) - lgamma(a) + a * theta - b * exp(theta);
}

/* gamma on the parameter x itself (INLA "gamma"): x ~ Gamma(shape = a, rate = b). */
static inline double gmrfs_prior_gamma(double x, double a, double b)
{
	return a * log(b) - lgamma(a) + (a - 1.0) * log(x) - b * x;
}

/* gaussian / normal on theta (INLA "gaussian"): param = (mean, prec).
 * prec == 0 is treated as an (improper) flat prior. */
static inline double gmrfs_prior_gaussian(double theta, double mean, double prec)
{
	if (prec == 0.0) return 0.0;
	double d = theta - mean;
	return -0.5 * log(2.0 * M_PI) + 0.5 * log(prec) - 0.5 * prec * d * d;
}

/* improper flat prior on theta (INLA "flat"). */
static inline double gmrfs_prior_flat(double theta)
{
	(void) theta;
	return 0.0;
}

/* ---- dispatchers -------------------------------------------------------- */

/* Integer prior codes (passed from R as an int -- robust across INLA
 * versions, unlike character arguments). */
#define GMRFS_PRIOR_HORSESHOE 1
#define GMRFS_PRIOR_PC_PREC   2
#define GMRFS_PRIOR_LOGGAMMA  3
#define GMRFS_PRIOR_GAMMA     4
#define GMRFS_PRIOR_GAUSSIAN  5
#define GMRFS_PRIOR_FLAT      6

/* dispatch by integer code; NAN if unknown / too few params. */
static inline double gmrfs_log_prior_code(int code, double theta,
					  const double *p, int np)
{
	switch (code) {
	case GMRFS_PRIOR_HORSESHOE: return (np >= 1) ? gmrfs_prior_horseshoe(theta, p[0]) : NAN;
	case GMRFS_PRIOR_PC_PREC:   return (np >= 2) ? gmrfs_prior_pc_prec(theta, p[0], p[1]) : NAN;
	case GMRFS_PRIOR_LOGGAMMA:  return (np >= 2) ? gmrfs_prior_loggamma(theta, p[0], p[1]) : NAN;
	case GMRFS_PRIOR_GAMMA:     return (np >= 2) ? gmrfs_prior_gamma(theta, p[0], p[1]) : NAN;
	case GMRFS_PRIOR_GAUSSIAN:  return (np >= 2) ? gmrfs_prior_gaussian(theta, p[0], p[1]) : NAN;
	case GMRFS_PRIOR_FLAT:      return gmrfs_prior_flat(theta);
	default:                    return NAN;
	}
}

/* Returns the log prior density for the named prior; NAN if the name is
 * unknown or too few parameters were supplied. */
static inline double gmrfs_log_prior(const char *name, double theta,
				     const double *p, int np)
{
	if (!name) return NAN;
	if (!strcasecmp(name, "horseshoe") || !strcasecmp(name, "hs"))
		return (np >= 1) ? gmrfs_prior_horseshoe(theta, p[0]) : NAN;
	if (!strcasecmp(name, "pc.prec") || !strcasecmp(name, "pcprec"))
		return (np >= 2) ? gmrfs_prior_pc_prec(theta, p[0], p[1]) : NAN;
	if (!strcasecmp(name, "loggamma"))
		return (np >= 2) ? gmrfs_prior_loggamma(theta, p[0], p[1]) : NAN;
	if (!strcasecmp(name, "gamma"))
		return (np >= 2) ? gmrfs_prior_gamma(theta, p[0], p[1]) : NAN;
	if (!strcasecmp(name, "gaussian") || !strcasecmp(name, "normal"))
		return (np >= 2) ? gmrfs_prior_gaussian(theta, p[0], p[1]) : NAN;
	if (!strcasecmp(name, "flat"))
		return gmrfs_prior_flat(theta);
	return NAN;
}

/* ---- shared block-model log prior --------------------------------------- */

/* Block-model log-prior helper shared by cblocks, pblocks, sblocks.
 *
 * Handles three orthogonal extensions on top of a per-block horseshoe:
 *   * noise_scale: read the Gaussian family's log-precision (appended
 *     to `theta` beyond the block's own p entries, as INLA does when
 *     data->theta_all_names lists it) and multiply
 *     sigma_e = exp(-theta_gauss/2) into the horseshoe scale, so
 *     sigma_j ~ C+(0, tau0 * sigma_e).
 *   * marg_tau: integrate the global scale tau out of the prior via a
 *     60-point log-tau quadrature with a half-Cauchy(0,1) hyperprior.
 *     tau never surfaces as an INLA hyperparameter.
 *   * plain (both above off): reproduce the pre-existing behaviour of
 *     the block models, i.e. gmrfs_log_prior_code with (tau0, ...) or
 *     the user-supplied param.
 *
 * Params:
 *   data          -- cgeneric data pointer (used to read `tau0`,
 *                    `noise_scale`, `marg_tau`).
 *   theta         -- theta vector from INLA. Entries 0..p-1 are the
 *                    block's own log-precisions. When noise_scale is
 *                    on, entries p..(p + data->theta_all_n - 1) carry
 *                    the "other" hyperparameters INLA is tracking,
 *                    in the order given by data->theta_all_names.
 *   p             -- number of blocks (length of the block's own
 *                    theta prefix).
 *   prior_code    -- prior id from the priors.h catalog. Used only in
 *                    the plain path (marg_tau = 0, noise_scale = 0).
 *   prior_params  -- user prior parameters; if NULL, defaults to (tau0).
 *                    Used only in the plain path.
 *   prior_nparams -- length of prior_params (if not NULL).
 *
 * Returns log p(theta_1..p [| sigma_e]).
 */
static inline double gmrfs_log_prior_blocks(inla_cgeneric_data_tp *data,
					    const double *theta, int p,
					    int prior_code,
					    const double *prior_params,
					    int prior_nparams)
{
	/* ---- read options --------------------------------------------- */
	double tau0 = 1.0;
	int noise_scale = 0, marg_tau = 0;
	for (int k = 0; k < data->n_doubles; k++) {
		if (!strcasecmp(data->doubles[k]->name, "tau0")) {
			tau0 = data->doubles[k]->doubles[0];
		}
	}
	for (int k = 0; k < data->n_ints; k++) {
		if (!strcasecmp(data->ints[k]->name, "noise_scale"))
			noise_scale = data->ints[k]->ints[0];
		else if (!strcasecmp(data->ints[k]->name, "marg_tau"))
			marg_tau = data->ints[k]->ints[0];
	}

	/* ---- read sigma_e from theta[p + i] (only if noise_scale) ----- */
	double sigma_e = 1.0;
	if (noise_scale && theta && data->theta_all_names) {
		for (int i = 0; i < data->theta_all_n; i++) {
			if (data->theta_all_names[i] &&
			    strstr(data->theta_all_names[i], "Gaussian")) {
				sigma_e = exp(-0.5 * theta[p + i]);
				break;
			}
		}
	}

	/* ---- marg_tau path: 60-point log-tau quadrature --------------- */
	if (marg_tau) {
		const int K = 60;
		const double g_lo = -4.0, g_hi = 4.0;
		double dg = (g_hi - g_lo) / (double) (K - 1);
		double t[60];
		double mx = -1e300;
		for (int k = 0; k < K; k++) {
			double g = g_lo + k * dg;
			double tau_k = exp(g);
			double te = noise_scale ? tau_k * sigma_e : tau_k;
			double Sk = 0.0;
			for (int j = 0; j < p; j++)
				Sk += gmrfs_prior_horseshoe(theta[j], te);
			/* log pi(tau_k) + Jacobian + log(dg):
			 *   log(2/pi) - log(1+tau^2)   +   g   +   log(dg) */
			t[k] = Sk +
			       (log(2.0 / M_PI) - log(1.0 + tau_k * tau_k)) +
			       g +
			       log(dg);
			if (t[k] > mx) mx = t[k];
		}
		double s = 0.0;
		for (int k = 0; k < K; k++) s += exp(t[k] - mx);
		return mx + log(s);
	}

	/* ---- plain path: sum of per-block priors ---------------------- */
	double tau_eff = noise_scale ? tau0 * sigma_e : tau0;
	double def_p[1] = { tau_eff };
	const double *pp; int np;
	if (noise_scale) {
		/* noise_scale + non-marg-tau forces the horseshoe form
		 * because the tau_eff replacement only makes sense for
		 * horseshoe. */
		pp = def_p; np = 1;
		prior_code = GMRFS_PRIOR_HORSESHOE;
	} else if (prior_params) {
		pp = prior_params; np = prior_nparams;
	} else {
		pp = def_p; np = 1;
	}
	double lp = 0.0;
	for (int j = 0; j < p; j++)
		lp += gmrfs_log_prior_code(prior_code, theta[j], pp, np);
	return lp;
}

#endif
