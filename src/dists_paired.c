#include "dists_paired.h"

//' R_haversine_paired
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_paired (SEXP x_, SEXP y_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    double *rx, *ry, *rout;
    double cosy1, cosy2;

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < n; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        cosy1 = cos (rx [n + i] * M_PI / 180.0); // y-value of x data
        cosy2 = cos (ry [n + i] * M_PI / 180.0);
        rout [i] = one_haversine (rx [i], rx [n + i], ry [i], ry [n + i],
                cosy1, cosy2);
    }

    UNPROTECT (3);

    return out;
}

//' R_vincenty_paired
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_paired (SEXP x_, SEXP y_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    double *rx, *ry, *rout;
    double siny1, cosy1, siny2, cosy2;

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < n; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        siny1 = sin (rx [n + i] * M_PI / 180.0); // y-value of x data
        cosy1 = cos (rx [n + i] * M_PI / 180.0); // y-value of x data
        siny2 = sin (ry [n + i] * M_PI / 180.0);
        cosy2 = cos (ry [n + i] * M_PI / 180.0);
        rout [i] = one_vincenty (rx [i], ry [i],
                siny1, cosy1, siny2, cosy2);
    }

    UNPROTECT (3);

    return out;
}

//' R_cheap_paired
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_paired (SEXP x_, SEXP y_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    double *rx, *ry, *rout;
    double ymin = 9999.9, ymax = -9999.9;
    double cosy;

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    // Get maximal latitude range
    for (size_t i = 0; i < n; i++)
    {
        if (rx [n + i] < ymin)
            ymin = rx [n + i];
        if (rx [n + i] > ymax)
            ymax = rx [n + i];
        if (ry [n + i] < ymin)
            ymin = ry [n + i];
        if (ry [n + i] > ymax)
            ymax = ry [n + i];
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    cosy = cos ((ymin + ymax) / 2.0);

    for (size_t i = 0; i < n; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        rout [i] = one_cheap (rx [i], rx [n + i], ry [i], ry [n + i], cosy);
    }

    UNPROTECT (3);

    return out;
}


//' R_geodesic_paired
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_paired (SEXP x_, SEXP y_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    double *rx, *ry, *rout;

    SEXP out = PROTECT (allocVector (REALSXP, n));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < n; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        rout [i] = one_geodesic (rx [i], rx [n + i], ry [i], ry [n + i]);
    }

    UNPROTECT (3);

    return out;
}
