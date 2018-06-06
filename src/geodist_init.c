#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _geodist_rcpp_haversine(SEXP);
extern SEXP _geodist_rcpp_haversine_xy(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_geodist_rcpp_haversine",    (DL_FUNC) &_geodist_rcpp_haversine,    1},
    {"_geodist_rcpp_haversine_xy", (DL_FUNC) &_geodist_rcpp_haversine_xy, 2},
    {NULL, NULL, 0}
};

void R_init_geodist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
