#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

//' R_haversine
//' @param x_ Single vector of x-values
//' @param x_ Single vector of y-values
//' @noRd
SEXP R_haversine_vec (SEXP x_, SEXP y_)
{
    size_t n = length (x_);
    size_t n2 = n * n;
    //Rprintf ("n = %d ; len = %d \n", n, n2);

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    double cosy1 [n]; // y-values are indexed in [n+1:n]
    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (ry [i] * M_PI / 180.0);
        rout [i * n + i] = 0.0;
    }

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = (i + 1); j < n; j++)
        {
            size_t indx1 = i * n + j;
            size_t indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_haversine (rx [i], ry [i],
                    rx [j], ry [j], cosy1 [i], cosy1 [j]);
        }
    }

    UNPROTECT (3);

    return out;
}

//' R_vincenty
//' @param x_ Single vector of x-values
//' @param x_ Single vector of y-values
//' @noRd
SEXP R_vincenty_vec (SEXP x_, SEXP y_)
{
    size_t n = length (x_);
    size_t n2 = n * n;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    double siny1 [n], cosy1 [n]; // y-values are indexed in [n+1:n]
    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (ry [i] * M_PI / 180.0);
        siny1 [i] = sin (ry [i] * M_PI / 180.0);
        rout [i * n + i] = 0.0;
    }

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = (i + 1); j < n; j++)
        {
            size_t indx1 = i * n + j;
            size_t indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_vincenty (rx [i], ry [i],
                    rx [j], ry [j],
                    siny1 [i], cosy1 [i], siny1 [j], cosy1 [j]);
        }
    }

    UNPROTECT (3);

    return out;
}

//' R_cheap
//' @param x_ Single vector of x-values
//' @param x_ Single vector of y-values
//' @noRd
SEXP R_cheap_vec (SEXP x_, SEXP y_)
{
    size_t n = length (x_);
    size_t n2 = n * n;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    // Get maximal latitude range
    double ymin = 9999.9, ymax = -9999.9;
    for (size_t i = 0; i < n; i++)
    {
        if (ry [i] < ymin)
            ymin = ry [i];
        if (ry [i] > ymax)
            ymax = ry [i];
        rout [i * n + i] = 0.0;
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    double cosy = cos ((ymin + ymax) / 2.0);

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = (i + 1); j < n; j++)
        {
            size_t indx1 = i * n + j;
            size_t indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_cheap (rx [i], ry [i],
                    rx [j], ry [j], cosy);
        }
    }

    UNPROTECT (3);

    return out;
}

//' R_geodesic
//' @param x_ Single vector of x-values
//' @param x_ Single vector of y-values
//' @noRd
SEXP R_geodesic_vec (SEXP x_, SEXP y_)
{
    size_t n = length (x_);
    size_t n2 = n * n;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < n; i++)
        rout [i * n + i] = 0.0;

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = (i + 1); j < n; j++)
        {
            size_t indx1 = i * n + j;
            size_t indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_geodesic (rx [i], ry [i],
                    rx [j], ry [j]);
        }
    }

    UNPROTECT (3);

    return out;
}
