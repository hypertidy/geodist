#include <R.h>
#include <Rinternals.h>

#include "common.h"
#include "WSG84-defs.h"

//' R_haversine_seq
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_seq (SEXP x_)
{
    size_t n = floor (length (x_) / 2);

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

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

    UNPROTECT (2);

    return out;
}

//' R_vincenty_seq
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_seq (SEXP x_)
{
    size_t n = floor (length (x_) / 2);

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

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

    UNPROTECT (2);

    return out;
}

//' R_cheap_seq
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_seq (SEXP x_)
{
    size_t n = floor (length (x_) / 2);

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    // Get maximal latitude range
    double ymin = 9999.9, ymax = -9999.9;
    for (size_t i = 0; i < n; i++)
    {
        if (rx [n + i] < ymin)
            ymin = rx [n + i];
        if (rx [n + i] > ymax)
            ymax = rx [n + i];
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    double cosy = cos ((ymin + ymax) / 2.0);

    rout [0] = NA_REAL;

    for (size_t i = 1; i < n; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        rout [i] = one_cheap (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i], cosy);
    }

    UNPROTECT (2);

    return out;
}

//' R_geodesic_seq
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_seq (SEXP x_)
{
    size_t n = floor (length (x_) / 2);

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    rout [0] = NA_REAL;

    for (size_t i = 1; i < n; i++)
    {
        rout [i] = one_geodesic (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i]);
    }

    UNPROTECT (2);

    return out;
}
