#include <R.h>
#include <Rinternals.h>

#include "common.h"
#include "WSG84-defs.h"

//' R_haversine_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_seq_range (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    double *rx;
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    rx = REAL (x_);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 1; i < n; i++)
    {
        double cosy1 = cos (rx [n + i] * M_PI / 180.0);
        double d = one_haversine (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i],
                cos (rx [n + i] * M_PI / 180.0),
                cos (rx [n + i - 1] * M_PI / 180.0));
        if (d < min)
            min = d;
        if (d > max)
            max = d;
    }

    double *rout;
    SEXP out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (2);

    return out;
}

//' R_vincenty_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_seq_range (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    double *rx;
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    rx = REAL (x_);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 1; i < n; i++)
    {
        double siny1 = sin (rx [n + i - 1] * M_PI / 180.0);
        double cosy1 = cos (rx [n + i - 1] * M_PI / 180.0);
        double siny2 = sin (rx [n + i] * M_PI / 180.0);
        double cosy2 = cos (rx [n + i] * M_PI / 180.0);
        double d = one_vincenty (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i], siny1, cosy1, siny2, cosy2);
        if (d < min)
            min = d;
        if (d > max)
            max = d;
    }

    double *rout;
    SEXP out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (2);

    return out;
}

//' R_cheap_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_seq_range (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    double *rx;
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    rx = REAL (x_);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 1; i < n; i++)
    {
        double cosy = cos ((rx [i - 1] * M_PI / 180.0 +
                            rx [i] * M_PI / 180.0) / 2.0);
        double d = one_cheap (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i], cosy);
        if (d < min)
            min = d;
        if (d > max)
            max = d;
    }

    double *rout;
    SEXP out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (2);

    return out;
}

//' R_geodesic_seq_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_seq_range (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    double *rx;
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    rx = REAL (x_);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 1; i < n; i++)
    {
        double d = one_geodesic (rx [i - 1], rx [n + i - 1],
                rx [i], rx [n + i]);
        if (d < min)
            min = d;
        if (d > max)
            max = d;
    }

    double *rout;
    SEXP out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (2);

    return out;
}
