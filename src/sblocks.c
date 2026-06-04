/*
 * sblocks.c -- block-horseshoe random-effect selection cgeneric (stacked).
 *
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"
#include "priors.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Calloc(n_, type_) (type_ *) calloc((n_), sizeof(type_))

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

/* prior parameter vector */
static double *find_doubles_named(inla_cgeneric_data_tp *data, const char *name, int *len)
{
    for (int k = 0; k < data->n_doubles; k++) {
        if (!strcasecmp(data->doubles[k]->name, name)) {
            if (len) *len = data->doubles[k]->len;
            return data->doubles[k]->doubles;
        }
    }
    if (len) *len = 0;
    return NULL;
}

/* effective per-block log-precisions: estimated `theta`, or fixed values
 * (recycled from length 1) when fixed != 0. Caller frees. */
static double *gmrfs_theta_eff(int fixed, const double *theta_fixed, int lf,
                               const double *theta, int p)
{
    double *th = Calloc(p, double);
    for (int j = 0; j < p; j++)
        th[j] = fixed ? theta_fixed[(lf == 1) ? 0 : j] : (theta ? theta[j] : 0.0);
    return th;
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
double *inla_cgeneric_gmrfs_sblocks(inla_cgeneric_cmd_tp cmd,
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

    /* fixed mode: hold block precisions at theta_fixed, expose 0 hyperpar */
    int fixed = find_int_named(data, "fixed");
    if (fixed < 0) fixed = 0;
    int lf = 0;
    double *theta_fixed = find_doubles_named(data, "theta_fixed", &lf);
    double *th_eff = gmrfs_theta_eff(fixed, theta_fixed, lf, theta, p);

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
            double w = exp(th_eff[j]);
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
        if (fixed) {
            ret = Calloc(1, double);
            ret[0] = 0.0;                    /* no hyperparameters -- precisions fixed */
        } else {
            ret = Calloc(1 + p, double);
            ret[0] = (double) p;
            for (int j = 0; j < p; j++) ret[1 + j] = 1.0;
        }
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
        for (int j = 0; j < p; j++) sum_theta += th_eff[j];
        double half_log_det_Q = 0.5 * ((double) block_n * sum_theta + logdet_R_all);
        ret = Calloc(1, double);
        ret[0] = half_log_det_Q - 0.5 * (double) N * log(2.0 * M_PI);
        break;
    }

    case INLA_CGENERIC_LOG_PRIOR: {
        ret = Calloc(1, double);
        if (fixed) {
            ret[0] = 0.0;                    /* no hyperparameters -> no prior */
        } else {
            /* prior from the catalog (priors.h), per block log-precision;
             * defaults to the horseshoe with scale tau0. */
            int code = find_int_named(data, "prior_id");
            if (code < 0) code = GMRFS_PRIOR_HORSESHOE;
            int np = 0;
            double *pp = find_doubles_named(data, "param", &np);
            double defp[1] = { find_double_named_default(data, "tau0", 1.0) };
            if (!pp) { pp = defp; np = 1; }
            double lp = 0.0;
            for (int j = 0; j < p; j++) lp += gmrfs_log_prior_code(code, theta[j], pp, np);
            ret[0] = lp;
        }
        break;
    }

    case INLA_CGENERIC_VOID:
    case INLA_CGENERIC_QUIT:
    default:
        break;
    }

    free(th_eff);
    return ret;
}
