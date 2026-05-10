/*
 * cblocks.c -- block-horseshoe random-effect selection cgeneric (collapsed).
 *
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "cgeneric.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Calloc(n_, type_) (type_ *) calloc((n_), sizeof(type_))

/* LAPACK */
extern void dpotrf_(const char *uplo, const int *n, double *a, const int *lda, int *info);
extern void dpotri_(const char *uplo, const int *n, double *a, const int *lda, int *info);

/* ------------------------------------------------------------------ */
/*  Numerically stable e^x * E_1(x), x > 0.                            */
/* ------------------------------------------------------------------ */
static double expE1(double x)
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
    static const double LOG_PI_SQRT_2PI = 1.8378770664093453;
    double tau02 = tau0 * tau0;
    double x = (0.5 / tau02) * exp(-theta);
    double v = expE1(x);
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

/* Build Sigma = sum_j exp(-theta_j) G_j into the upper triangle of dst
   (column-major, n x n).  Only the upper is read by dpotrf("U", ...). */
static void assemble_sigma_upper(double *dst, int n, int p,
                                 const double *G_rowmajor,
                                 const double *theta)
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
}

static double log_det_from_upper_chol(const double *Sigma, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++)
        s += 2.0 * log(Sigma[i + (size_t) i * n]);
    return s;
}

/* ================================================================== */
double *inla_cgeneric_cblocks(inla_cgeneric_cmd_tp cmd,
                                   double *theta,
                                   inla_cgeneric_data_tp *data)
{
    double *ret = NULL;

    assert(!strcasecmp(data->ints[0]->name, "n"));
    int n = data->ints[0]->ints[0];          /* latent dim = block size */
    assert(n > 0);

    int p = find_int_named(data, "p");
    assert(p > 0);

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
        assemble_sigma_upper(Sigma, n, p, G->x, theta);

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
        ret = Calloc(1 + p, double);
        ret[0] = (double) p;
        for (int j = 0; j < p; j++) ret[1 + j] = 1.0;
        break;
    }

    case INLA_CGENERIC_LOG_NORM_CONST: {
        /* c(theta) = -0.5 log|Sigma(theta)| - 0.5 n log(2 pi).
           Stateless: factor a fresh Sigma scratch every call. */
        double *Sigma = Calloc((size_t) n * n, double);
        assemble_sigma_upper(Sigma, n, p, G->x, theta);
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
