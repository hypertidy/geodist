#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

//' R_haversine_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_range (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    double *rx;
    rx = REAL (x_);

    double cosy1 [n]; // y-values are indexed in [n+1:n]
    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (rx [n + i] * M_PI / 180.0);
    }

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = (i + 1); j < n; j++)
        {
            double d = one_haversine (rx [i], rx [n + i],
                    rx [j], rx [n + j], cosy1 [i], cosy1 [j]);
            if (d < min)
                min = d;
            else if (d > max)
                max = d;
        }
    }

    double *rout;
    SEXP out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (1);

    return out;
}

//' R_vincenty_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_range (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    size_t n2 = n * n;
    double *rx;
    rx = REAL (x_);

    double siny1 [n], cosy1 [n]; // y-values are indexed in [n+1:n]
    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (rx [n + i] * M_PI / 180.0);
        siny1 [i] = sin (rx [n + i] * M_PI / 180.0);
    }

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = (i + 1); j < n; j++)
        {
            double d = one_vincenty (rx [i], rx [n + i],
                    rx [j], rx [n + j],
                    siny1 [i], cosy1 [i], siny1 [j], cosy1 [j]);
            if (d < min)
                min = d;
            else if (d > max)
                max = d;
        }
    }

    double *rout;
    SEXP out = PROTECT (allocVector (REALSXP, n2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (1);

    return out;
}

//' R_cheap_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_range (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    double *rx;
    rx = REAL (x_);

    // Get maximal latitude range
    double ymin = 9999.9, ymax = -9999.9;
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
    double cosy = cos ((ymin + ymax) / 2.0);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = (i + 1); j < n; j++)
        {
            double d = one_cheap (rx [i], rx [n + i],
                    rx [j], rx [n + j], cosy);
            if (d < min)
                min = d;
            else if (d > max)
                max = d;
        }
    }

    double *rout;
    SEXP out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (1);

    return out;
}

//' R_geodesic_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_range (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    double *rx;
    rx = REAL (x_);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = (i + 1); j < n; j++)
        {
            double d = one_geodesic (rx [i], rx [n + i],
                    rx [j], rx [n + j]);
            if (d < min)
                min = d;
            else if (d > max)
                max = d;
        }
    }

    double *rout;
    SEXP out = PROTECT (allocVector (REALSXP, 2));
    rout = REAL (out);
    rout [0] = min;
    rout [1] = max;

    UNPROTECT (1);

    return out;
}

