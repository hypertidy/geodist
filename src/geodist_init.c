#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP R_cheap(SEXP);
extern SEXP R_cheap_paired(SEXP, SEXP);
extern SEXP R_cheap_paired_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_cheap_range(SEXP);
extern SEXP R_cheap_seq(SEXP);
extern SEXP R_cheap_seq_range(SEXP);
extern SEXP R_cheap_seq_vec(SEXP, SEXP);
extern SEXP R_cheap_vec(SEXP, SEXP);
extern SEXP R_cheap_xy(SEXP, SEXP);
extern SEXP R_cheap_xy_min(SEXP, SEXP);
extern SEXP R_cheap_xy_range(SEXP, SEXP);
extern SEXP R_cheap_xy_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_geodesic(SEXP);
extern SEXP R_geodesic_paired(SEXP, SEXP);
extern SEXP R_geodesic_paired_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_geodesic_range(SEXP);
extern SEXP R_geodesic_seq(SEXP);
extern SEXP R_geodesic_seq_range(SEXP);
extern SEXP R_geodesic_seq_vec(SEXP, SEXP);
extern SEXP R_geodesic_vec(SEXP, SEXP);
extern SEXP R_geodesic_xy(SEXP, SEXP);
extern SEXP R_geodesic_xy_min(SEXP, SEXP);
extern SEXP R_geodesic_xy_range(SEXP, SEXP);
extern SEXP R_geodesic_xy_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_haversine(SEXP);
extern SEXP R_haversine_paired(SEXP, SEXP);
extern SEXP R_haversine_paired_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_haversine_range(SEXP);
extern SEXP R_haversine_seq(SEXP);
extern SEXP R_haversine_seq_range(SEXP);
extern SEXP R_haversine_seq_vec(SEXP, SEXP);
extern SEXP R_haversine_vec(SEXP, SEXP);
extern SEXP R_haversine_xy(SEXP, SEXP);
extern SEXP R_haversine_xy_min(SEXP, SEXP);
extern SEXP R_haversine_xy_range(SEXP, SEXP);
extern SEXP R_haversine_xy_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_vincenty(SEXP);
extern SEXP R_vincenty_paired(SEXP, SEXP);
extern SEXP R_vincenty_paired_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_vincenty_range(SEXP);
extern SEXP R_vincenty_seq(SEXP);
extern SEXP R_vincenty_seq_range(SEXP);
extern SEXP R_vincenty_seq_vec(SEXP, SEXP);
extern SEXP R_vincenty_vec(SEXP, SEXP);
extern SEXP R_vincenty_xy(SEXP, SEXP);
extern SEXP R_vincenty_xy_min(SEXP, SEXP);
extern SEXP R_vincenty_xy_range(SEXP, SEXP);
extern SEXP R_vincenty_xy_vec(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_cheap",                (DL_FUNC) &R_cheap,                1},
    {"R_cheap_paired",         (DL_FUNC) &R_cheap_paired,         2},
    {"R_cheap_paired_vec",     (DL_FUNC) &R_cheap_paired_vec,     4},
    {"R_cheap_range",          (DL_FUNC) &R_cheap_range,          1},
    {"R_cheap_seq",            (DL_FUNC) &R_cheap_seq,            1},
    {"R_cheap_seq_range",      (DL_FUNC) &R_cheap_seq_range,      1},
    {"R_cheap_seq_vec",        (DL_FUNC) &R_cheap_seq_vec,        2},
    {"R_cheap_vec",            (DL_FUNC) &R_cheap_vec,            2},
    {"R_cheap_xy",             (DL_FUNC) &R_cheap_xy,             2},
    {"R_cheap_xy_min",         (DL_FUNC) &R_cheap_xy_min,         2},
    {"R_cheap_xy_range",       (DL_FUNC) &R_cheap_xy_range,       2},
    {"R_cheap_xy_vec",         (DL_FUNC) &R_cheap_xy_vec,         4},
    {"R_geodesic",             (DL_FUNC) &R_geodesic,             1},
    {"R_geodesic_paired",      (DL_FUNC) &R_geodesic_paired,      2},
    {"R_geodesic_paired_vec",  (DL_FUNC) &R_geodesic_paired_vec,  4},
    {"R_geodesic_range",       (DL_FUNC) &R_geodesic_range,       1},
    {"R_geodesic_seq",         (DL_FUNC) &R_geodesic_seq,         1},
    {"R_geodesic_seq_range",   (DL_FUNC) &R_geodesic_seq_range,   1},
    {"R_geodesic_seq_vec",     (DL_FUNC) &R_geodesic_seq_vec,     2},
    {"R_geodesic_vec",         (DL_FUNC) &R_geodesic_vec,         2},
    {"R_geodesic_xy",          (DL_FUNC) &R_geodesic_xy,          2},
    {"R_geodesic_xy_min",      (DL_FUNC) &R_geodesic_xy_min,      2},
    {"R_geodesic_xy_range",    (DL_FUNC) &R_geodesic_xy_range,    2},
    {"R_geodesic_xy_vec",      (DL_FUNC) &R_geodesic_xy_vec,      4},
    {"R_haversine",            (DL_FUNC) &R_haversine,            1},
    {"R_haversine_paired",     (DL_FUNC) &R_haversine_paired,     2},
    {"R_haversine_paired_vec", (DL_FUNC) &R_haversine_paired_vec, 4},
    {"R_haversine_range",      (DL_FUNC) &R_haversine_range,      1},
    {"R_haversine_seq",        (DL_FUNC) &R_haversine_seq,        1},
    {"R_haversine_seq_range",  (DL_FUNC) &R_haversine_seq_range,  1},
    {"R_haversine_seq_vec",    (DL_FUNC) &R_haversine_seq_vec,    2},
    {"R_haversine_vec",        (DL_FUNC) &R_haversine_vec,        2},
    {"R_haversine_xy",         (DL_FUNC) &R_haversine_xy,         2},
    {"R_haversine_xy_min",     (DL_FUNC) &R_haversine_xy_min,     2},
    {"R_haversine_xy_range",   (DL_FUNC) &R_haversine_xy_range,   2},
    {"R_haversine_xy_vec",     (DL_FUNC) &R_haversine_xy_vec,     4},
    {"R_vincenty",             (DL_FUNC) &R_vincenty,             1},
    {"R_vincenty_paired",      (DL_FUNC) &R_vincenty_paired,      2},
    {"R_vincenty_paired_vec",  (DL_FUNC) &R_vincenty_paired_vec,  4},
    {"R_vincenty_range",       (DL_FUNC) &R_vincenty_range,       1},
    {"R_vincenty_seq",         (DL_FUNC) &R_vincenty_seq,         1},
    {"R_vincenty_seq_range",   (DL_FUNC) &R_vincenty_seq_range,   1},
    {"R_vincenty_seq_vec",     (DL_FUNC) &R_vincenty_seq_vec,     2},
    {"R_vincenty_vec",         (DL_FUNC) &R_vincenty_vec,         2},
    {"R_vincenty_xy",          (DL_FUNC) &R_vincenty_xy,          2},
    {"R_vincenty_xy_min",      (DL_FUNC) &R_vincenty_xy_min,      2},
    {"R_vincenty_xy_range",    (DL_FUNC) &R_vincenty_xy_range,    2},
    {"R_vincenty_xy_vec",      (DL_FUNC) &R_vincenty_xy_vec,      4},
    {NULL, NULL, 0}
};

void R_init_geodist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
