/*
 * qinv.c -- marginal variances diag(Q^-1) of a sparse precision matrix,
 * by the Takahashi recursions on its Cholesky factor (Rue and Martino 2007,
 * as R-INLA's inla.qinv()). Used to scale intrinsic models without INLA.
 *
 * Input: the simplicial LL' factor of Matrix::Cholesky(Q, LDL = FALSE,
 * super = FALSE): column j has nz[j] entries from p[j], the diagonal first.
 * The factor keeps its symbolic pattern, which is closed under the
 * recursions: if L_ki and L_ji are in the pattern (k > j > i), so is L_kj.
 * Output: diag(Sigma) in the order of the factor (permuted).
 */

/* R package only: INLA compiles every .c file in src/ without R, and with
 * GMRFS_R_PACKAGE unset (it is set in Makevars) this file is empty. */
#ifdef GMRFS_R_PACKAGE

#include <stdlib.h>

#include <R.h>
#include <Rinternals.h>

typedef struct {
    int row;
    double val;
} entry_tp;

static int cmp_entry(const void *a, const void *b)
{
    int ra = ((const entry_tp *) a)->row, rb = ((const entry_tp *) b)->row;
    return (ra > rb) - (ra < rb);
}

/* position of row r among the sorted off-diagonal rows of column c, or -1 */
static int find_row(const entry_tp *col, int len, int r)
{
    int lo = 0, hi = len - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        if (col[mid].row == r) return mid;
        if (col[mid].row < r) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

SEXP gmrfs_takahashi_diag(SEXP Sp, SEXP Si, SEXP Sx, SEXP Snz)
{
    int n = length(Snz);
    const int *p = INTEGER(Sp), *ii = INTEGER(Si), *nz = INTEGER(Snz);
    const double *x = REAL(Sx);

    /* off-diagonal entries of each column, sorted by row; Sigma on the same
     * pattern (S) and its diagonal (Sd) */
    int *start = (int *) R_alloc(n + 1, sizeof(int));
    double *Ld = (double *) R_alloc(n, sizeof(double));
    start[0] = 0;
    for (int j = 0; j < n; j++) {
        if (nz[j] < 1 || ii[p[j]] != j) error("the factor must have the diagonal first in each column");
        start[j + 1] = start[j] + nz[j] - 1;
        Ld[j] = x[p[j]];
    }
    entry_tp *L = (entry_tp *) R_alloc(start[n] > 0 ? start[n] : 1, sizeof(entry_tp));
    double *S = (double *) R_alloc(start[n] > 0 ? start[n] : 1, sizeof(double));
    for (int j = 0; j < n; j++) {
        for (int k = 1; k < nz[j]; k++) {
            L[start[j] + k - 1].row = ii[p[j] + k];
            L[start[j] + k - 1].val = x[p[j] + k];
        }
        qsort(L + start[j], nz[j] - 1, sizeof(entry_tp), cmp_entry);
    }

    SEXP ans = PROTECT(allocVector(REALSXP, n));
    double *Sd = REAL(ans);

    /* Sigma_kj, k, j >= the current column, from what is already computed */
#define SIGMA(k_, j_, out_)                                                   \
    do {                                                                      \
        int r_ = ((k_) > (j_) ? (k_) : (j_)), c_ = ((k_) > (j_) ? (j_) : (k_)); \
        if (r_ == c_) {                                                       \
            (out_) = Sd[c_];                                                  \
        } else {                                                              \
            int pos_ = find_row(L + start[c_], start[c_ + 1] - start[c_], r_); \
            if (pos_ < 0) error("the factor's pattern is not closed (entry %d, %d)", r_, c_); \
            (out_) = S[start[c_] + pos_];                                     \
        }                                                                     \
    } while (0)

    for (int j = n - 1; j >= 0; j--) {
        int m = start[j + 1] - start[j];
        const entry_tp *col = L + start[j];
        double *Scol = S + start[j];
        /* Sigma_{r_a, j} = -1/L_jj sum_b L_{r_b, j} Sigma_{r_b, r_a} */
        for (int a = m - 1; a >= 0; a--) {
            double sum = 0.0, s;
            for (int b = 0; b < m; b++) {
                SIGMA(col[b].row, col[a].row, s);
                sum += col[b].val * s;
            }
            Scol[a] = -sum / Ld[j];
        }
        /* Sigma_jj = 1/L_jj^2 - 1/L_jj sum_b L_{r_b, j} Sigma_{r_b, j} */
        double sum = 0.0;
        for (int b = 0; b < m; b++) sum += col[b].val * Scol[b];
        Sd[j] = 1.0 / (Ld[j] * Ld[j]) - sum / Ld[j];
    }
#undef SIGMA

    UNPROTECT(1);
    return ans;
}

#endif /* GMRFS_R_PACKAGE */
