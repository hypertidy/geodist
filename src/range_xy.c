#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

//' R_haversine_xy_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_xy_range (SEXP x_, SEXP y_)
{
    size_t nx = floor (length (x_) / 2);
    size_t ny = floor (length (y_) / 2);
    double *rx, *ry;
    rx = REAL (x_);
    ry = REAL (y_);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        double cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < ny; j++)
        {
            double cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            double d = one_haversine (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], cosy1, cosy2);
            if (d < min)
                min = d;
            if (d > max)
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

//' R_vincenty_xy_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_xy_range (SEXP x_, SEXP y_)
{
    size_t nx = floor (length (x_) / 2);
    size_t ny = floor (length (y_) / 2);
    double *rx, *ry;
    rx = REAL (x_);
    ry = REAL (y_);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        double siny1 = sin (rx [nx + i] * M_PI / 180.0); // y-value of x data
        double cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < ny; j++)
        {
            double siny2 = sin (ry [ny + j] * M_PI / 180.0);
            double cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            double d = one_vincenty (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], siny1, cosy1, siny2, cosy2);
            if (d < min)
                min = d;
            if (d > max)
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

//' R_cheap_xy_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_xy_range (SEXP x_, SEXP y_)
{
    size_t nx = floor (length (x_) / 2);
    size_t ny = floor (length (y_) / 2);
    double *rx, *ry;
    rx = REAL (x_);
    ry = REAL (y_);

    // Get maximal latitude range
    double ymin = 9999.9, ymax = -9999.9;
    for (size_t i = 0; i < nx; i++)
    {
        if (rx [nx + i] < ymin)
            ymin = rx [nx + i];
        if (rx [nx + i] > ymax)
            ymax = rx [nx + i];
    }
    for (size_t i = 0; i < ny; i++)
    {
        if (ry [ny + i] < ymin)
            ymin = ry [ny + i];
        if (ry [ny + i] > ymax)
            ymax = ry [ny + i];
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    double cosy = cos ((ymin + ymax) / 2.0);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = 0; j < ny; j++)
        {
            double d = one_cheap (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], cosy);
            if (d < min)
                min = d;
            if (d > max)
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


//' R_geodesic_xy_range
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_xy_range (SEXP x_, SEXP y_)
{
    size_t nx = floor (length (x_) / 2);
    size_t ny = floor (length (y_) / 2);
    double *rx, *ry;
    rx = REAL (x_);
    ry = REAL (y_);

    double min = 100.0 * equator, max = -100.0 * equator;

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 100 == 0)
            R_CheckUserInterrupt ();
        for (size_t j = 0; j < ny; j++)
        {
            double d = one_geodesic (rx [i], rx [nx + i],
                    ry [j], ry [ny + j]);
            if (d < min)
                min = d;
            if (d > max)
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
