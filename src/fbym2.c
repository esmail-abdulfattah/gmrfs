/*
 * fbym2.c -- flexible BYM2 cgeneric: bym2 whose structured part is fbesag.
 *
 *   v = (1/sqrt(tau)) * ( sqrt(1-phi) * w + sqrt(phi) * u )
 *   w ~ N(0, I),   u ~ N(0, Q_fb(d)^-)
 *
 * Q_fb(d) is the fbesag precision with relative partition weights
 * gamma_k = exp(d_k), sum_k d_k = 0: edge i~j has weight
 * c_i * (gamma_p(i) + gamma_p(j)) / 2, with c_i the scaling constant of the
 * connected component of i; an island has Q_fb_ii = gamma_p(i). With d = 0
 * this is INLA's bym2.
 *
 * Latent x = (v, u), length 2n. As in INLA's Qfunc_bym2:
 *   Q_vv = tau/(1-phi) I,  Q_vu = -sqrt(phi*tau)/(1-phi) I,
 *   Q_uu = phi/(1-phi) I + Q_fb(d).
 * A jitter EPS is added to the whole diagonal, so the null eigenvalues of Q
 * are exactly EPS for every theta and the log-determinant INLA computes stays
 * correct up to a constant. EPS also shifts each small non-null eigenvalue
 * lambda by log(1 + EPS/lambda), so it is kept small.
 *
 * theta = (log tau, logit phi, d_1, ..., d_{P-1}),  d_P = -(d_1 + ... + d_{P-1}).
 * Priors: pc.prec on tau; the bym2 PC prior on phi (tabulated in R on the
 * logit scale); d ~ N(0, sd^2 (I - J/P)) on the sum-to-zero subspace.
 *
 * The log prior also carries a constant K (from R) that removes the graph-only
 * terms INLA's bym2 leaves out of its normalising constant, so the mlik
 * follows bym2's convention. It does not change the posterior.
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"
#include "priors.h"

#define Calloc(n_, type_) (type_ *) calloc((n_), sizeof(type_))
#if 0
#define EPS 1e-5
#endif
/* EPS also shifts the small non-null eigenvalues; with 1e-5 the mlik on a
 * 150 x 150 lattice was 0.03 too high, with 1e-8 it is exact to 1e-3 */
#define EPS 1e-8

static inla_cgeneric_vec_tp *find_ivec(inla_cgeneric_data_tp *data, const char *name)
{
    for (int k = 0; k < data->n_ints; k++)
        if (!strcasecmp(data->ints[k]->name, name)) return data->ints[k];
    assert(0 && "missing integer argument");
    return NULL;
}

static inla_cgeneric_vec_tp *find_dvec(inla_cgeneric_data_tp *data, const char *name)
{
    for (int k = 0; k < data->n_doubles; k++)
        if (!strcasecmp(data->doubles[k]->name, name)) return data->doubles[k];
    assert(0 && "missing double argument");
    return NULL;
}

/* linear interpolation in a table with increasing x, linear extrapolation
 * from the end segments */
static double table_interp(const double *x, const double *y, int m, double x0)
{
    int lo, hi;
    if (x0 <= x[0]) {
        lo = 0; hi = 1;
    } else if (x0 >= x[m - 1]) {
        lo = m - 2; hi = m - 1;
    } else {
        lo = 0; hi = m - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x[mid] <= x0) lo = mid; else hi = mid;
        }
    }
    double t = (x0 - x[lo]) / (x[hi] - x[lo]);
    return y[lo] + t * (y[hi] - y[lo]);
}

/* log prior of logit(phi) beyond the table: the exact tail of the PC prior.
 * tail = (theta_max, d_max, r, lambda, c0): d^2 = d_max^2 + r (L(t) - L(theta_max)),
 * L(t) = log(1 + e^t), log prior = c0 - lambda d + log(r phi / d). */
static double phi_tail_log_prior(const double *tail, double t)
{
    double L = t + log1p(exp(-t));
    double Lmax = tail[0] + log1p(exp(-tail[0]));
    double d = sqrt(tail[1] * tail[1] + tail[2] * (L - Lmax));
    double phi = 1.0 / (1.0 + exp(-t));
    return tail[4] - tail[3] * d + log(tail[2] * phi / d);
}

double *inla_cgeneric_gmrfs_fbym2(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
    double *ret = NULL;

    int N = find_ivec(data, "n")->ints[0];             /* 2n */
    int n = N / 2;
    int npart = find_ivec(data, "npart")->ints[0];
    int ntheta = 2 + npart - 1;
    int *gr = find_ivec(data, "graph")->ints;           /* [len, N, M, ii..., jj...] */
    int M = gr[2];
    int *misc = find_ivec(data, "misc")->ints;
    double *scale = find_dvec(data, "scale")->doubles;

    double tau = 1.0, phi = 0.5, one_m_phi = 0.5, *gam = NULL;
    if (theta) {
        tau = exp(theta[0]);
        phi = 1.0 / (1.0 + exp(-theta[1]));
        one_m_phi = 1.0 / (1.0 + exp(theta[1]));   /* 1 - phi without cancellation */
        gam = Calloc(npart, double);
        double s = 0.0;
        for (int k = 0; k < npart - 1; k++) {
            gam[k] = exp(theta[2 + k]);
            s += theta[2 + k];
        }
        gam[npart - 1] = exp(-s);
    }

    switch (cmd) {

    case INLA_CGENERIC_VOID:
        assert(!(cmd == INLA_CGENERIC_VOID));
        break;

    case INLA_CGENERIC_GRAPH: {
        ret = Calloc(2 + 2 * M, double);
        ret[0] = N;
        ret[1] = M;
        for (int k = 0; k < M; k++) {
            ret[2 + k] = gr[3 + k];
            ret[2 + M + k] = gr[3 + M + k];
        }
        break;
    }

    case INLA_CGENERIC_Q: {
        /* same order as the graph: rows 0..n-1 give (i,i), (i,i+n); rows
         * n..2n-1 give (n+i,n+i) then (n+i,n+j) for neighbours j > i */
        if (0) {
            /* old: 1 - phi is 0 in double precision for logit(phi) > ~37 */
            double q_vv = tau / (1.0 - phi);
            double q_vu = -sqrt(phi * tau) / (1.0 - phi);
            double ridge = phi / (1.0 - phi);
            (void) q_vv; (void) q_vu; (void) ridge;
        }
        double q_vv = tau / one_m_phi;
        double q_vu = -sqrt(phi * tau) / one_m_phi;
        double ridge = phi / one_m_phi;

        ret = Calloc(2 + M, double);
        ret[0] = -1;
        ret[1] = M;
        int count = 2, k = 0;
        for (int i = 0; i < n; i++) {
            ret[count++] = q_vv + EPS;
            ret[count++] = q_vu;
        }
        for (int i = 0; i < n; i++) {
            int num_nei_i = misc[k++];
            double g_i = gam[misc[k++]];
            int size_neighbors = misc[k++];

            if (scale[i] < 0.0) {
                /* island */
                assert(num_nei_i == 0 && size_neighbors == 0);
                (void) num_nei_i;       /* used only by the assert */
                ret[count++] = ridge + g_i + EPS;
                continue;
            }

            int save_count = count++;
            double diag = 0.0;
            for (int j = 0; j < size_neighbors; j++) {
                double g_j = gam[misc[k++]];
                int tick = misc[k++];
                double w = scale[i] * 0.5 * (g_i + g_j);
                diag += w;
                if (tick > i) ret[count++] = -w;
            }
            ret[save_count] = ridge + diag + EPS;
        }
        assert(count == 2 + M);
        break;
    }

    case INLA_CGENERIC_MU: {
        ret = Calloc(1, double);
        ret[0] = 0;
        break;
    }

    case INLA_CGENERIC_INITIAL: {
        double *init = find_dvec(data, "initial")->doubles;
        ret = Calloc(ntheta + 1, double);
        ret[0] = ntheta;
        for (int i = 0; i < ntheta; i++) ret[1 + i] = init[i];
        break;
    }

    case INLA_CGENERIC_LOG_NORM_CONST:
        /* log|Q_fb(d)| has no closed form: INLA computes it */
        break;

    case INLA_CGENERIC_LOG_PRIOR: {
        double *prec = find_dvec(data, "prec")->doubles;           /* (u, alpha) */
        inla_cgeneric_vec_tp *px = find_dvec(data, "phi_x");
        double *py = find_dvec(data, "phi_y")->doubles;
        double sd = find_dvec(data, "sd_gamma")->doubles[0];

        double val = find_dvec(data, "K")->doubles[0];
        val += gmrfs_prior_pc_prec(theta[0], prec[0], prec[1]);
        if (0) {
            /* old: linear extrapolation on both sides; the right side rises,
             * so the prior was not proper */
            val += table_interp(px->doubles, py, px->len, theta[1]);
        }
        /* new: the left side is linear (slope ~1, the exact asymptote); the
         * right side has the exact tail */
        double *tail = find_dvec(data, "phi_tail")->doubles;
        if (theta[1] > tail[0]) {
            val += phi_tail_log_prior(tail, theta[1]);
        } else {
            val += table_interp(px->doubles, py, px->len, theta[1]);
        }

        if (npart > 1) {
            /* d on the sum-to-zero subspace, written in d_1..d_{P-1} */
            double ss = 0.0, s = 0.0;
            for (int k = 0; k < npart - 1; k++) {
                ss += theta[2 + k] * theta[2 + k];
                s += theta[2 + k];
            }
            ss += s * s;
            val += -0.5 * (npart - 1) * log(2.0 * M_PI * sd * sd) + 0.5 * log((double) npart)
                - 0.5 * ss / (sd * sd);
        }
        ret = Calloc(1, double);
        ret[0] = val;
        break;
    }

    case INLA_CGENERIC_QUIT:
    default:
        break;
    }

    free(gam);
    return ret;
}
