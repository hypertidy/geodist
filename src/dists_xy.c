#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

//' R_haversine_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_xy (SEXP x_, SEXP y_)
{
    size_t nx = floor (length (x_) / 2);
    size_t ny = floor (length (y_) / 2);
    size_t n2 = nx * ny;
    //Rprintf ("(nx, ny) = (%d , %d )\n", nx, ny);

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        double cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < ny; j++)
        {
            double cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            rout [i * ny + j] = one_haversine (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], cosy1, cosy2);
        }
    }

    UNPROTECT (3);

    return out;
}

//' R_vincenty_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_xy (SEXP x_, SEXP y_)
{
    size_t nx = floor (length (x_) / 2);
    size_t ny = floor (length (y_) / 2);
    size_t n2 = nx * ny;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        double siny1 = sin (rx [nx + i] * M_PI / 180.0); // y-value of x data
        double cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < ny; j++)
        {
            double siny2 = sin (ry [ny + j] * M_PI / 180.0);
            double cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            rout [i * ny + j] = one_vincenty (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], siny1, cosy1, siny2, cosy2);
        }
    }

    UNPROTECT (3);

    return out;
}

//' R_cheap_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_xy (SEXP x_, SEXP y_)
{
    size_t nx = floor (length (x_) / 2);
    size_t ny = floor (length (y_) / 2);
    size_t n2 = nx * ny;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

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

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = 0; j < ny; j++)
        {
            rout [i * ny + j] = one_cheap (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], cosy);
        }
    }

    UNPROTECT (3);

    return out;
}


//' R_geodesic_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_xy (SEXP x_, SEXP y_)
{
    size_t nx = floor (length (x_) / 2);
    size_t ny = floor (length (y_) / 2);
    size_t n2 = nx * ny;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = 0; j < ny; j++)
        {
            rout [i * ny + j] = one_geodesic (rx [i], rx [nx + i],
                    ry [j], ry [ny + j]);
        }
    }

    UNPROTECT (3);

    return out;
}
