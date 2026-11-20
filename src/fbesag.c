#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdbool.h>
#include "cgeneric.h"

#define Calloc(n_, type_) (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))

/* INLA's copy of this model is inla_cgeneric_fbesag, built into the inla
 * binary; a built-in model wins over any shlib with the same name, so the
 * package uses its own name. */
double *inla_cgeneric_gmrfs_fbesag(inla_cgeneric_cmd_tp cmd,
                                   double *theta,
                                   inla_cgeneric_data_tp *data)
{
    double *ret = NULL, *prec = NULL;

    assert(!strcasecmp(data->ints[0]->name, "n"));
    int N = data->ints[0]->ints[0];
    assert(N > 0);

    int M = data->ints[3]->ints[2];
    int npart = data->ints[2]->ints[0];

    if (theta) {
        prec = Calloc(npart, double);
        for (int i = 0; i < npart; i++) prec[i] = exp(theta[i]);
    }

    switch (cmd) {

    case INLA_CGENERIC_VOID:
        assert(!(cmd == INLA_CGENERIC_VOID));
        break;

    case INLA_CGENERIC_GRAPH: {
        ret = Calloc(2 + 2 * M, double);
        ret[0] = data->ints[3]->ints[1];   /* graph dim N */
        ret[1] = M;                         /* nnz */
        for (int i = 0; i < M; i++) {
            ret[2 + i]     = data->ints[3]->ints[3 + i];
            ret[2 + M + i] = data->ints[3]->ints[3 + M + i];
        }
        break;
    }

    case INLA_CGENERIC_Q: {
        double scaled_cnst = data->doubles[1]->doubles[2];
        ret = Calloc(2 + M, double);
        ret[0] = -1;
        ret[1] = M;
        int count = 2, k = 0;
        if (0) {
        /* old: one scaling constant for the whole graph */
        for (int i = 0; i < N; i++) {
            int num_nei_i = data->ints[4]->ints[k]; k++;
            double mas_prec = prec[data->ints[4]->ints[k]]; k++;
            int size_neighbors = data->ints[4]->ints[k]; k++;

            int save_count = count;
            ret[count] = 0;
            count++;
            double tmp_neighbors = 0.0, sum_tmp_neighbors = 0.0;
            for (int j = 0; j < size_neighbors; j++) {
                tmp_neighbors = -0.5 * prec[data->ints[4]->ints[k]]; k++;
                int tick = data->ints[4]->ints[k]; k++;
                if (tick > i) {
                    ret[count] = scaled_cnst * (tmp_neighbors - 0.5 * mas_prec);
                    count++;
                }
                sum_tmp_neighbors += tmp_neighbors;
            }
            ret[save_count] = scaled_cnst * (-sum_tmp_neighbors + 0.5 * num_nei_i * mas_prec) + 1e-5;
        }
        }

        /* new: as Qfunc_besag() in INLA. scale[i] is the constant of the
         * connected component of node i; an edge never leaves its
         * component. scale[i] < 0 marks an island: Q_ii = tau_p(i). */
        assert(!strcasecmp(data->doubles[3]->name, "scale"));
        double *scale = data->doubles[3]->doubles;
        for (int i = 0; i < N; i++) {
            int num_nei_i = data->ints[4]->ints[k]; k++;
            double mas_prec = prec[data->ints[4]->ints[k]]; k++;
            int size_neighbors = data->ints[4]->ints[k]; k++;

            if (scale[i] < 0.0) {
                assert(num_nei_i == 0 && size_neighbors == 0);
                ret[count] = mas_prec;
                count++;
                continue;
            }

            int save_count = count;
            ret[count] = 0;
            count++;
            double tmp_neighbors = 0.0, sum_tmp_neighbors = 0.0;
            for (int j = 0; j < size_neighbors; j++) {
                tmp_neighbors = -0.5 * prec[data->ints[4]->ints[k]]; k++;
                int tick = data->ints[4]->ints[k]; k++;
                if (tick > i) {
                    ret[count] = scale[i] * (tmp_neighbors - 0.5 * mas_prec);
                    count++;
                }
                sum_tmp_neighbors += tmp_neighbors;
            }
            if (0) ret[save_count] = scale[i] * (-sum_tmp_neighbors + 0.5 * num_nei_i * mas_prec) + 1e-5;
            /* a small jitter: it also shifts the small non-null eigenvalues,
             * which with 1e-5 biased the mlik on large graphs */
            ret[save_count] = scale[i] * (-sum_tmp_neighbors + 0.5 * num_nei_i * mas_prec) + 1e-8;
        }
        break;
    }

    case INLA_CGENERIC_MU: {
        ret = Calloc(npart, double);
        for (int i = 0; i < npart; i++) ret[i] = 0;
        break;
    }

    case INLA_CGENERIC_INITIAL: {
        ret = Calloc(npart + 1, double);
        ret[0] = npart;
        for (int i = 1; i <= npart; i++)
            ret[i] = data->doubles[2]->doubles[i - 1];
        break;
    }

    case INLA_CGENERIC_LOG_NORM_CONST:
        /* ret stays NULL: INLA computes it (break, so prec is freed) */
        break;

    case INLA_CGENERIC_LOG_PRIOR: {
        double lam = data->doubles[0]->doubles[0], mean_theta = 0;
        double val1 = data->doubles[1]->doubles[0], val2 = data->doubles[1]->doubles[1];

        for (int i = 0; i < npart; i++) mean_theta += theta[i];
        mean_theta = mean_theta / npart;

        double res1 = log(lam) - lam * exp(-0.5 * mean_theta) - 0.5 * mean_theta;
        double res2 = 0;
        for (int i = 0; i < npart; i++) {
            res2 += val1 * SQR(theta[i] - mean_theta);
            for (int j = 0; j < i; j++)
                res2 += 2 * val2 * (theta[i] - mean_theta) * (theta[j] - mean_theta);
        }

        ret = Calloc(1, double);
        if (0) {
        /* old: correct up to a constant only */
        ret[0] = res1 - 0.5 * res2;
        }

        /* new: a proper density in theta, so the mlik can be compared with
         * other models. theta = mean * 1 + x, x on the sum-to-zero subspace:
         * pc.prec on the mean (with its 1/2), N(0, sd^2) on x within the
         * subspace, and 1/sqrt(P) from the change of coordinates. K (from R)
         * puts the normalising constant of the field on besag's convention. */
        assert(!strcasecmp(data->doubles[4]->name, "sd_gamma"));
        assert(!strcasecmp(data->doubles[5]->name, "K"));
        double sd = data->doubles[4]->doubles[0];
        double val = res1 + log(0.5) + data->doubles[5]->doubles[0];
        if (npart > 1)
            val += -0.5 * (npart - 1) * log(2.0 * M_PI * sd * sd) - 0.5 * log((double) npart) - 0.5 * res2;
        ret[0] = val;
        break;
    }

    case INLA_CGENERIC_QUIT:
    default:
        break;
    }

    free(prec);
    return ret;
}
