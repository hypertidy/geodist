#include "range_x.h"

//' R_haversine_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_range (SEXP x_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    double *rx;

    double cosy1 [n]; // y-values are indexed in [n+1:n]
    double min = 100.0 * equator, max = -100.0 * equator;
    double d, *rout;
    SEXP out;

    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    rx = REAL (x_);

    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (rx [n + i] * M_PI / 180.0);
    }

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = (i + 1); j < n; j++)
        {
            d = one_haversine (rx [i], rx [n + i],
                    rx [j], rx [n + j], cosy1 [i], cosy1 [j]);
            if (d < min)
                min = d;
            if (d > max)
                max = d;
        }
    }

    out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (2);

    return out;
}

//' R_vincenty_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_range (SEXP x_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    size_t n2 = n * n;
    double *rx;

    double siny1 [n], cosy1 [n]; // y-values are indexed in [n+1:n]
    double min = 100.0 * equator, max = -100.0 * equator;
    double d, *rout;
    SEXP out;

    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    rx = REAL (x_);

    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (rx [n + i] * M_PI / 180.0);
        siny1 [i] = sin (rx [n + i] * M_PI / 180.0);
    }

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = (i + 1); j < n; j++)
        {
            d = one_vincenty (rx [i], rx [j],
                    siny1 [i], cosy1 [i], siny1 [j], cosy1 [j]);
            if (d < min)
                min = d;
            if (d > max)
                max = d;
        }
    }

    out = PROTECT (allocVector (REALSXP, n2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (2);

    return out;
}

//' R_cheap_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_range (SEXP x_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    double *rx;

    double ymin = 9999.9, ymax = -9999.9;
    double min = 100.0 * equator, max = -100.0 * equator;
    double cosy, d, *rout;
    SEXP out;

    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    rx = REAL (x_);

    // Get maximal latitude range
    for (size_t i = 0; i < n; i++)
    {
        if (rx [n + i] < ymin)
            ymin = rx [n + i];
        else if (rx [n + i] > ymax)
            ymax = rx [n + i];
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    cosy = cos ((ymin + ymax) / 2.0);

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = (i + 1); j < n; j++)
        {
            d = one_cheap (rx [i], rx [n + i],
                    rx [j], rx [n + j], cosy);
            if (d < min)
                min = d;
            if (d > max)
                max = d;
        }
    }

    out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (2);

    return out;
}

//' R_geodesic_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_range (SEXP x_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    double *rx;

    double min = 100.0 * equator, max = -100.0 * equator;
    double d, *rout;
    SEXP out;

    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    rx = REAL (x_);

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = (i + 1); j < n; j++)
        {
            d = one_geodesic (rx [i], rx [n + i],
                    rx [j], rx [n + j]);
            if (d < min)
                min = d;
            if (d > max)
                max = d;
        }
    }

    out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (2);

    return out;
}
