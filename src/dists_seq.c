#include <R.h>
#include <Rinternals.h>

#include "common.h"
#include "WSG84-defs.h"

//' R_haversine
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_seq (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    SEXP out = PROTECT (allocVector (REALSXP, n));
    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    rout [0] = NA_REAL;
    for (size_t i = 1; i < n; i++)
    {
        double cosy1 = cos (rx [n + i] * M_PI / 180.0);
        rout [i] = one_haversine (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i],
                cos (rx [n + i] * M_PI / 180.0),
                cos (rx [n + i - 1] * M_PI / 180.0));
    }

    UNPROTECT (1);

    return out;
}

//' R_vincenty
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_seq (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    SEXP out = PROTECT (allocVector (REALSXP, n));
    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    rout [0] = NA_REAL;

    for (size_t i = 1; i < n; i++)
    {
        double siny1 = sin (rx [n + i - 1] * M_PI / 180.0);
        double cosy1 = cos (rx [n + i - 1] * M_PI / 180.0);
        double siny2 = sin (rx [n + i] * M_PI / 180.0);
        double cosy2 = cos (rx [n + i] * M_PI / 180.0);
        rout [i] = one_vincenty (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i], siny1, cosy1, siny2, cosy2);
    }

    UNPROTECT (1);

    return out;
}

//' R_cheap
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_seq (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    SEXP out = PROTECT (allocVector (REALSXP, n));
    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    rout [0] = NA_REAL;

    for (size_t i = 1; i < n; i++)
    {
        double cosy = cos ((rx [i - 1] * M_PI / 180.0 +
                            rx [i] * M_PI / 180.0) / 2.0);
        rout [i] = one_cheap (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i], cosy);
    }

    UNPROTECT (1);

    return out;
}

//' R_geodesic_seq
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_seq (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    SEXP out = PROTECT (allocVector (REALSXP, n));
    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    rout [0] = NA_REAL;

    for (size_t i = 1; i < n; i++)
    {
        rout [i] = one_geodesic (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i]);
    }

    UNPROTECT (1);

    return out;
}

