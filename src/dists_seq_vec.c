#include "dists_seq_vec.h"

//' R_haversine_seq_vec
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_seq_vec (SEXP x_, SEXP y_)
{
    size_t n = (size_t) length (x_);
    double *rx, *ry, *rout;

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    rout [0] = NA_REAL;
    for (size_t i = 1; i < n; i++)
    {
        rout [i] = one_haversine (rx [i - 1], ry [i - 1],
                rx [i], ry [i],
                cos (ry [i] * M_PI / 180.0),
                cos (ry [i - 1] * M_PI / 180.0));
    }

    UNPROTECT (2);

    return out;
}

//' R_vincenty_seq_vec
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_seq_vec (SEXP x_, SEXP y_)
{
    size_t n = (size_t) length (x_);
    double *rx, *ry, *rout;
    double siny1, cosy1, siny2, cosy2;

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    rout [0] = NA_REAL;

    for (size_t i = 1; i < n; i++)
    {
        siny1 = sin (ry [i - 1] * M_PI / 180.0);
        cosy1 = cos (ry [i - 1] * M_PI / 180.0);
        siny2 = sin (ry [i] * M_PI / 180.0);
        cosy2 = cos (ry [i] * M_PI / 180.0);
        rout [i] = one_vincenty (rx [i - 1], rx [i],
                siny1, cosy1, siny2, cosy2);
    }

    UNPROTECT (2);

    return out;
}

//' R_cheap_seq_vec
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_seq_vec (SEXP x_, SEXP y_)
{
    size_t n = (size_t) length (x_);
    double *rx, *ry, *rout;
    double ymin = 9999.9, ymax = -9999.9;
    double cosy;

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    // Get maximal latitude range
    for (size_t i = 0; i < n; i++)
    {
        if (ry [i] < ymin)
            ymin = ry [i];
        if (ry [i] > ymax)
            ymax = ry [i];
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    cosy = cos ((ymin + ymax) / 2.0);

    rout [0] = NA_REAL;

    for (size_t i = 1; i < n; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        rout [i] = one_cheap (rx [i - 1], ry [i - 1],
                rx [i], ry [i], cosy);
    }

    UNPROTECT (2);

    return out;
}

//' R_geodesic_seq_vec
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_seq_vec (SEXP x_, SEXP y_)
{
    size_t n = (size_t) length (x_);
    double *rx, *ry, *rout;

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    rout [0] = NA_REAL;

    for (size_t i = 1; i < n; i++)
    {
        rout [i] = one_geodesic (rx [i - 1], ry [i - 1],
                rx [i], ry [i]);
    }

    UNPROTECT (2);

    return out;
}
