/*
 * init.c -- R-level entry point for the cgeneric models.
 *
 * gmrfs_cgeneric_eval() builds an inla_cgeneric_data_tp from the R list
 * that inla.cgeneric.define() stores (ints, doubles, characters, matrices,
 * smatrices) and calls the model function directly, so the models can be
 * evaluated from R without INLA.
 */

/* R package only: INLA compiles every .c file in src/ without R, and with
 * GMRFS_R_PACKAGE unset (it is set in Makevars) this file is empty. */
#ifdef GMRFS_R_PACKAGE

#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "cgeneric.h"

inla_cgeneric_func_tp inla_cgeneric_gmrfs_fbesag;
inla_cgeneric_func_tp inla_cgeneric_gmrfs_fbym2;
inla_cgeneric_func_tp inla_cgeneric_gmrfs_cblocks;
inla_cgeneric_func_tp inla_cgeneric_gmrfs_pblocks;
inla_cgeneric_func_tp inla_cgeneric_gmrfs_sblocks;

static inla_cgeneric_func_tp *find_model(const char *name)
{
    if (!strcmp(name, "inla_cgeneric_gmrfs_fbesag")) return inla_cgeneric_gmrfs_fbesag;
    if (!strcmp(name, "inla_cgeneric_gmrfs_fbym2")) return inla_cgeneric_gmrfs_fbym2;
    if (!strcmp(name, "inla_cgeneric_gmrfs_cblocks")) return inla_cgeneric_gmrfs_cblocks;
    if (!strcmp(name, "inla_cgeneric_gmrfs_pblocks")) return inla_cgeneric_gmrfs_pblocks;
    if (!strcmp(name, "inla_cgeneric_gmrfs_sblocks")) return inla_cgeneric_gmrfs_sblocks;
    return NULL;
}

static char *elt_name(SEXP names, int k)
{
    const char *s = CHAR(STRING_ELT(names, k));
    char *out = R_alloc(strlen(s) + 1, sizeof(char));
    strcpy(out, s);
    return out;
}

/* all storage is R_alloc'ed, so it is released when .Call returns */
static inla_cgeneric_data_tp *build_data(SEXP data)
{
    SEXP ints = VECTOR_ELT(data, 0), doubles = VECTOR_ELT(data, 1), chars = VECTOR_ELT(data, 2),
        mats = VECTOR_ELT(data, 3), smats = VECTOR_ELT(data, 4);

    inla_cgeneric_data_tp *d = (inla_cgeneric_data_tp *) R_alloc(1, sizeof(inla_cgeneric_data_tp));
    memset(d, 0, sizeof(inla_cgeneric_data_tp));
    d->threads.max = d->threads.outer = d->threads.inner = 1;

    d->n_ints = LENGTH(ints);
    d->ints = (inla_cgeneric_vec_tp **) R_alloc(d->n_ints + 1, sizeof(inla_cgeneric_vec_tp *));
    for (int k = 0; k < d->n_ints; k++) {
        SEXP x = VECTOR_ELT(ints, k);
        inla_cgeneric_vec_tp *v = (inla_cgeneric_vec_tp *) R_alloc(1, sizeof(inla_cgeneric_vec_tp));
        memset(v, 0, sizeof(inla_cgeneric_vec_tp));
        v->name = elt_name(getAttrib(ints, R_NamesSymbol), k);
        v->len = LENGTH(x);
        v->ints = INTEGER(x);
        d->ints[k] = v;
    }

    d->n_doubles = LENGTH(doubles);
    d->doubles = (inla_cgeneric_vec_tp **) R_alloc(d->n_doubles + 1, sizeof(inla_cgeneric_vec_tp *));
    for (int k = 0; k < d->n_doubles; k++) {
        SEXP x = VECTOR_ELT(doubles, k);
        inla_cgeneric_vec_tp *v = (inla_cgeneric_vec_tp *) R_alloc(1, sizeof(inla_cgeneric_vec_tp));
        memset(v, 0, sizeof(inla_cgeneric_vec_tp));
        v->name = elt_name(getAttrib(doubles, R_NamesSymbol), k);
        v->len = LENGTH(x);
        v->doubles = REAL(x);
        d->doubles[k] = v;
    }

    d->n_chars = LENGTH(chars);
    d->chars = (inla_cgeneric_vec_tp **) R_alloc(d->n_chars + 1, sizeof(inla_cgeneric_vec_tp *));
    for (int k = 0; k < d->n_chars; k++) {
        inla_cgeneric_vec_tp *v = (inla_cgeneric_vec_tp *) R_alloc(1, sizeof(inla_cgeneric_vec_tp));
        memset(v, 0, sizeof(inla_cgeneric_vec_tp));
        v->name = elt_name(getAttrib(chars, R_NamesSymbol), k);
        v->chars = elt_name(VECTOR_ELT(chars, k), 0);
        v->len = (int) strlen(v->chars);
        d->chars[k] = v;
    }

    /* dense: c(nrow, ncol, x), x stored row by row */
    d->n_mats = LENGTH(mats);
    d->mats = (inla_cgeneric_mat_tp **) R_alloc(d->n_mats + 1, sizeof(inla_cgeneric_mat_tp *));
    for (int k = 0; k < d->n_mats; k++) {
        double *x = REAL(VECTOR_ELT(mats, k));
        inla_cgeneric_mat_tp *m = (inla_cgeneric_mat_tp *) R_alloc(1, sizeof(inla_cgeneric_mat_tp));
        m->name = elt_name(getAttrib(mats, R_NamesSymbol), k);
        m->nrow = (int) x[0];
        m->ncol = (int) x[1];
        m->x = x + 2;
        d->mats[k] = m;
    }

    /* sparse: c(nrow, ncol, n, i, j, x), 0-based triplets */
    d->n_smats = LENGTH(smats);
    d->smats = (inla_cgeneric_smat_tp **) R_alloc(d->n_smats + 1, sizeof(inla_cgeneric_smat_tp *));
    for (int k = 0; k < d->n_smats; k++) {
        double *x = REAL(VECTOR_ELT(smats, k));
        inla_cgeneric_smat_tp *m = (inla_cgeneric_smat_tp *) R_alloc(1, sizeof(inla_cgeneric_smat_tp));
        m->name = elt_name(getAttrib(smats, R_NamesSymbol), k);
        m->nrow = (int) x[0];
        m->ncol = (int) x[1];
        m->n = (int) x[2];
        m->i = (int *) R_alloc(m->n + 1, sizeof(int));
        m->j = (int *) R_alloc(m->n + 1, sizeof(int));
        for (int t = 0; t < m->n; t++) {
            m->i[t] = (int) x[3 + t];
            m->j[t] = (int) x[3 + m->n + t];
        }
        m->x = x + 3 + 2 * m->n;
        d->smats[k] = m;
    }

    return d;
}

/* length of the array a cgeneric model returns for `cmd` */
static R_xlen_t ret_length(inla_cgeneric_cmd_tp cmd, const double *ret)
{
    switch (cmd) {
    case INLA_CGENERIC_GRAPH:
        return 2 + 2 * (R_xlen_t) ret[1];
    case INLA_CGENERIC_Q:
        /* optimised format (-1, M, x) or the full one (n, M, i, j, x) */
        return (ret[0] == -1.0 ? 2 + (R_xlen_t) ret[1] : 2 + 3 * (R_xlen_t) ret[1]);
    case INLA_CGENERIC_MU:
    case INLA_CGENERIC_INITIAL:
        return 1 + (R_xlen_t) ret[0];
    case INLA_CGENERIC_LOG_NORM_CONST:
    case INLA_CGENERIC_LOG_PRIOR:
        return 1;
    default:
        return 0;
    }
}

/* .Call entry: returns the raw array of the model, or NULL if it returns none */
SEXP gmrfs_cgeneric_eval(SEXP model, SEXP cmd, SEXP theta, SEXP data)
{
    inla_cgeneric_func_tp *fun = find_model(CHAR(STRING_ELT(model, 0)));
    if (!fun) error("unknown cgeneric model '%s'", CHAR(STRING_ELT(model, 0)));

    inla_cgeneric_cmd_tp c = (inla_cgeneric_cmd_tp) asInteger(cmd);
    if (c < INLA_CGENERIC_Q || c > INLA_CGENERIC_LOG_PRIOR) error("invalid cgeneric command");

    inla_cgeneric_data_tp *d = build_data(data);
    double *ret = fun(c, (LENGTH(theta) > 0 ? REAL(theta) : NULL), d);
    if (!ret) return R_NilValue;

    R_xlen_t len = ret_length(c, ret);
    SEXP ans = PROTECT(allocVector(REALSXP, len));
    if (len > 0) memcpy(REAL(ans), ret, (size_t) len * sizeof(double));
    free(ret);
    UNPROTECT(1);
    return ans;
}

/* qinv.c */
SEXP gmrfs_takahashi_diag(SEXP Sp, SEXP Si, SEXP Sx, SEXP Snz);

static const R_CallMethodDef call_methods[] = {
    {"gmrfs_cgeneric_eval", (DL_FUNC) &gmrfs_cgeneric_eval, 4},
    {"gmrfs_takahashi_diag", (DL_FUNC) &gmrfs_takahashi_diag, 4},
    {NULL, NULL, 0}
};

void R_init_gmrfs(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

#endif /* GMRFS_R_PACKAGE */
