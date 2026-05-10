/*
 * stacked.c
 *
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Calloc(n_, type_) (type_ *) calloc((n_), sizeof(type_))

/* ------------------------------------------------------------------ */
/*  Numerically stable e^x * E_1(x), x > 0.                           */
/*  Small x : convergent series.                                       */
/*  Large x : even continued fraction via modified Lentz.              */
/* ------------------------------------------------------------------ */
static double expE1(double x)
{
    if (!(x > 0.0)) return NAN;

    if (x < 1.0) {
        /* E_1(x) = -gamma - log x + sum_{k>=1} (-1)^{k+1} x^k / (k*k!) */
        const double EG = 0.5772156649015329;
        double e1 = -EG - log(x);
        double term = 1.0;            /* term = (-x)^k / k! */
        for (int k = 1; k < 100; k++) {
            term *= -x / (double) k;
            double add = -term / (double) k;
            e1 += add;
            if (fabs(add) < 1e-16 * (1.0 + fabs(e1))) break;
        }
        return exp(x) * e1;
    } else {
        /* e^x E_1(x) = 1/(x+1 - 1^2/(x+3 - 2^2/(x+5 - ...))) */
        const double TINY = 1e-30;
        double b = x + 1.0;
        double d = 1.0 / b;
        double c = 1.0 / TINY;
        double h = d;
        for (int i = 1; i < 200; i++) {
            double a = -(double)(i * i);
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

static double log_prior_theta(double theta, double tau0)
{
    static const double LOG_PI_SQRT_2PI = 1.8378770664093453; /* log(pi*sqrt(2*pi)) */
    double tau02 = tau0 * tau0;
    double x = (0.5 / tau02) * exp(-theta);
    double v = expE1(x);
    /* v is strictly positive for x > 0; safe to log */
    return -0.5 * theta + log(v) - LOG_PI_SQRT_2PI - 2.0 * log(tau0);
}

/* ------------------------------------------------------------------ */
/*  Lookup helpers                                                     */
/* ------------------------------------------------------------------ */
static int find_int_named(inla_cgeneric_data_tp *data, const char *name)
{
    for (int k = 0; k < data->n_ints; k++) {
        if (!strcasecmp(data->ints[k]->name, name)) {
            return data->ints[k]->ints[0];
        }
    }
    return -1;
}

static double find_double_named_default(inla_cgeneric_data_tp *data,
                                        const char *name, double dflt)
{
    for (int k = 0; k < data->n_doubles; k++) {
        if (!strcasecmp(data->doubles[k]->name, name)) {
            return data->doubles[k]->doubles[0];
        }
    }
    return dflt;
}

static inla_cgeneric_mat_tp *find_mat_named(inla_cgeneric_data_tp *data,
                                            const char *name)
{
    for (int k = 0; k < data->n_mats; k++) {
        if (!strcasecmp(data->mats[k]->name, name)) {
            return data->mats[k];
        }
    }
    return NULL;
}

/* ================================================================== */
/* stacked: latent dim p*n, block-diagonal Q with diagonal blocks
 *   exp(theta_j) * R_j,   R_j a dense n x n positive-definite matrix.
 *
 * The blocks are passed dense as `R_stacked`, an (p*n) x n row-major
 * matrix (block j at rows j*n..(j+1)*n-1).  Graph + Q are emitted in
 * row-major upper-triangular order per block -- the SAME pattern that
 * partial.c uses for its bottom-right block, which INLA handles
 * correctly.  This is the load-bearing structural choice: emitting in
 * column-major (the order inla.as.sparse produces) makes INLA's inner
 * Newton-Raphson on the latent diverge.
 */
double *inla_cgeneric_stacked(inla_cgeneric_cmd_tp cmd,
                               double *theta,
                               inla_cgeneric_data_tp *data)
{
    double *ret = NULL;

    /* ---- mandatory args ------------------------------------------- */
    assert(!strcasecmp(data->ints[0]->name, "n"));
    int N = data->ints[0]->ints[0];
    assert(N > 0);

    int p = find_int_named(data, "p");
    assert(p > 0);
    assert(N % p == 0);
    int block_n = N / p;

    inla_cgeneric_mat_tp *R_stacked = find_mat_named(data, "R_stacked");
    assert(R_stacked != NULL);
    assert(R_stacked->nrow == p * block_n);
    assert(R_stacked->ncol == block_n);
    /* Block j at offset j*block_n*block_n in R_stacked->x (row-major). */

    int M_per_block = block_n * (block_n + 1) / 2;
    int M           = p * M_per_block;

    switch (cmd) {

    case INLA_CGENERIC_GRAPH: {
        ret = Calloc(2 + 2 * M, double);
        ret[0] = (double) N;
        ret[1] = (double) M;
        int k = 0;
        for (int j = 0; j < p; j++) {
            int offset = j * block_n;
            for (int a = 0; a < block_n; a++)
                for (int b = a; b < block_n; b++) {
                    ret[2 + k]     = (double) (offset + a);
                    ret[2 + M + k] = (double) (offset + b);
                    k++;
                }
        }
        break;
    }

    case INLA_CGENERIC_Q: {
        ret = Calloc(2 + M, double);
        ret[0] = -1.0;           /* optimised format (same order as graph) */
        ret[1] = (double) M;
        int k = 0;
        for (int j = 0; j < p; j++) {
            double w = exp(theta[j]);
            const double *Rj = R_stacked->x + (size_t) j * block_n * block_n;
            /* Rj is row-major n x n: entry (a, b) at Rj[a*block_n + b]. */
            for (int a = 0; a < block_n; a++)
                for (int b = a; b < block_n; b++) {
                    ret[2 + k] = w * Rj[(size_t) a * block_n + b];
                    k++;
                }
        }
        break;
    }

    case INLA_CGENERIC_MU: {
        ret = Calloc(1, double);
        ret[0] = 0.0;            /* zero-mean */
        break;
    }

    case INLA_CGENERIC_INITIAL: {
        ret = Calloc(1 + p, double);
        ret[0] = (double) p;
        for (int j = 0; j < p; j++) ret[1 + j] = 1.0;
        break;
    }

    case INLA_CGENERIC_LOG_NORM_CONST: {
        /* Closed form when logdet_R_all = sum_j log|R_j| is supplied
           from R: log|Q| = block_n * sum theta_j + logdet_R_all.
           Falls back to INLA's own factorisation if not supplied
           (which is fragile when R_j is ill-conditioned). */
        double logdet_R_all = find_double_named_default(data, "logdet_R_all", NAN);
        if (isnan(logdet_R_all)) {
            ret = NULL;
            break;
        }
        double sum_theta = 0.0;
        for (int j = 0; j < p; j++) sum_theta += theta[j];
        double half_log_det_Q = 0.5 * ((double) block_n * sum_theta + logdet_R_all);
        ret = Calloc(1, double);
        ret[0] = half_log_det_Q - 0.5 * (double) N * log(2.0 * M_PI);
        break;
    }

    case INLA_CGENERIC_LOG_PRIOR: {
        double tau0 = find_double_named_default(data, "tau0", 1.0);
        double lp = 0.0;
        for (int j = 0; j < p; j++) lp += log_prior_theta(theta[j], tau0);
        ret = Calloc(1, double);
        ret[0] = lp;
        break;
    }

    case INLA_CGENERIC_VOID:
    case INLA_CGENERIC_QUIT:
    default:
        break;
    }

    return ret;
}
