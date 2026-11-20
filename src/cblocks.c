/*
 * cblocks.c -- block-horseshoe random-effect selection cgeneric (collapsed).
 *
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "cgeneric.h"
#include "priors.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Calloc(n_, type_) (type_ *) calloc((n_), sizeof(type_))

/* LAPACK */
extern void dpotrf_(const char *uplo, const int *n, double *a, const int *lda, int *info);
extern void dpotri_(const char *uplo, const int *n, double *a, const int *lda, int *info);

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

/* effective per-block log-precisions: the estimated `theta`, or the supplied
 * fixed values (recycled from length 1) when fixed != 0. Caller frees. */
static double *gmrfs_theta_eff(int fixed, const double *theta_fixed, int lf,
                               const double *theta, int p)
{
    double *th = Calloc(p, double);
    for (int j = 0; j < p; j++)
        th[j] = fixed ? theta_fixed[(lf == 1) ? 0 : j] : (theta ? theta[j] : 0.0);
    return th;
}

/* Build Sigma = sum_j exp(-theta_j) G_j + nugget * I into the upper
   triangle of dst (column-major, n x n).  Only the upper is read by
   dpotrf("U", ...).  The optional `nugget` term is a diagonal floor
   that keeps Sigma well-conditioned when all blocks are strongly
   shrunk; pass 0.0 to reproduce the pre-nugget behaviour. */
static void assemble_sigma_upper(double *dst, int n, int p,
                                 const double *G_rowmajor,
                                 const double *theta, double nugget)
{
    for (int b = 0; b < n; b++)
        for (int a = 0; a <= b; a++)
            dst[a + (size_t) b * n] = 0.0;

    for (int j = 0; j < p; j++) {
        double w = exp(-theta[j]);
        const double *Gj = G_rowmajor + (size_t) j * n * n;
        for (int a = 0; a < n; a++)
            for (int b = a; b < n; b++)
                dst[a + (size_t) b * n] += w * Gj[(size_t) a * n + b];
    }
    if (nugget > 0.0) {
        for (int i = 0; i < n; i++)
            dst[i + (size_t) i * n] += nugget;
    }
}

static double log_det_from_upper_chol(const double *Sigma, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++)
        s += 2.0 * log(Sigma[i + (size_t) i * n]);
    return s;
}

/* ================================================================== */
double *inla_cgeneric_gmrfs_cblocks(inla_cgeneric_cmd_tp cmd,
                                   double *theta,
                                   inla_cgeneric_data_tp *data)
{
    double *ret = NULL;

    assert(!strcasecmp(data->ints[0]->name, "n"));
    int n = data->ints[0]->ints[0];          /* latent dim = block size */
    assert(n > 0);

    int p = find_int_named(data, "p");
    assert(p > 0);
    if (n <= 0 || p <= 0) return NULL;

    /* fixed mode: hold the block precisions at theta_fixed, expose 0 hyperpar */
    int fixed = find_int_named(data, "fixed");
    if (fixed < 0) fixed = 0;
    int lf = 0;
    double *theta_fixed = find_doubles_named(data, "theta_fixed", &lf);

    inla_cgeneric_mat_tp *G = find_mat_named(data, "G");
    assert(G != NULL);
    assert(G->nrow == n * p);
    assert(G->ncol == n);
    /* G->x is row-major: block j occupies rows j*n .. (j+1)*n-1,
       G_{j+1}[a,b] = G->x[(j*n + a) * n + b]. */

    int M = n * (n + 1) / 2;                 /* upper-triangle size     */

    switch (cmd) {

    case INLA_CGENERIC_GRAPH: {
        ret = Calloc(2 + 2 * M, double);
        ret[0] = (double) n;
        ret[1] = (double) M;
        int k = 0;
        for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++) {
                ret[2 + k]     = (double) i;
                ret[2 + M + k] = (double) j;
                k++;
            }
        break;
    }

    case INLA_CGENERIC_Q: {
        /* Stateless: allocate Sigma scratch per call. Trivially thread-safe. */
        double *Sigma = Calloc((size_t) n * n, double);
        double *th = gmrfs_theta_eff(fixed, theta_fixed, lf, theta, p);
        double nugget = find_double_named_default(data, "nugget", 0.0);
        assemble_sigma_upper(Sigma, n, p, G->x, th, nugget);
        free(th);

        int info = 0;
        dpotrf_("U", &n, Sigma, &n, &info);
        if (info == 0) dpotri_("U", &n, Sigma, &n, &info);
        if (info != 0) {
            free(Sigma);
            ret = Calloc(2 + M, double);
            ret[0] = -1.0;
            ret[1] = (double) M;
            for (int k = 0; k < M; k++) ret[2 + k] = NAN;
            break;
        }

        ret = Calloc(2 + M, double);
        ret[0] = -1.0;
        ret[1] = (double) M;
        int k = 0;
        for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++) {
                ret[2 + k] = Sigma[i + (size_t) j * n];
                k++;
            }
        free(Sigma);
        break;
    }

    case INLA_CGENERIC_MU: {
        ret = Calloc(1, double);
        ret[0] = 0.0;
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
        /* c(theta) = -0.5 log|Sigma(theta)| - 0.5 n log(2 pi).
           Stateless: factor a fresh Sigma scratch every call. */
        double *Sigma = Calloc((size_t) n * n, double);
        double *th = gmrfs_theta_eff(fixed, theta_fixed, lf, theta, p);
        double nugget = find_double_named_default(data, "nugget", 0.0);
        assemble_sigma_upper(Sigma, n, p, G->x, th, nugget);
        free(th);
        int info = 0;
        dpotrf_("U", &n, Sigma, &n, &info);
        if (info != 0) {
            free(Sigma);
            ret = Calloc(1, double);
            ret[0] = NAN;
            break;
        }
        double log_det = log_det_from_upper_chol(Sigma, n);
        free(Sigma);
        ret = Calloc(1, double);
        ret[0] = -0.5 * log_det - 0.5 * (double) n * log(2.0 * M_PI);
        break;
    }

    case INLA_CGENERIC_LOG_PRIOR: {
        ret = Calloc(1, double);
        if (fixed) {
            ret[0] = 0.0;                    /* no hyperparameters -> no prior */
        } else {
            /* Delegate to the shared block-model helper (priors.h):
             *   - plain path: sum of per-block priors from the catalog
             *   - noise_scale: multiplies sigma_e into the horseshoe scale
             *   - marg_tau:    60-point log-tau quadrature, tau integrated out
             * All three routes read their options from `data` directly. */
            int code = find_int_named(data, "prior_id");
            if (code < 0) code = GMRFS_PRIOR_HORSESHOE;
            int np = 0;
            double *pp = find_doubles_named(data, "param", &np);
            ret[0] = gmrfs_log_prior_blocks(data, theta, p, code, pp, np);
        }
        break;
    }

    case INLA_CGENERIC_VOID:
    case INLA_CGENERIC_QUIT:
    default:
        break;
    }

    return ret;
}
