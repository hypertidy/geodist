#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

//' R_haversine_xy_vec
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_xy_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_)
{
    size_t n1 = length (x1_);
    size_t n2 = length (x2_);
    size_t nm = n1 * n2;
    SEXP out = PROTECT (allocVector (REALSXP, nm));
    double *rx1, *ry1, *rx2, *ry2, *rout;
    rx1 = REAL (x1_);
    ry1 = REAL (y1_);
    rx2 = REAL (x2_);
    ry2 = REAL (y2_);
    rout = REAL (out);

    for (size_t i = 0; i < n1; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        double cosy1 = cos (ry1 [i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < n2; j++)
        {
            double cosy2 = cos (ry2 [j] * M_PI / 180.0);
            rout [i * n2 + j] = one_haversine (rx1 [i], ry1 [i],
                    rx2 [j], ry2 [j], cosy1, cosy2);
        }
    }

    UNPROTECT (1);

    return out;
}

//' R_vincenty_xy_vec
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_xy_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_)
{
    size_t n1 = length (x1_);
    size_t n2 = length (x2_);
    size_t nm = n1 * n2;
    SEXP out = PROTECT (allocVector (REALSXP, nm));
    double *rx1, *ry1, *rx2, *ry2, *rout;
    rx1 = REAL (x1_);
    ry1 = REAL (y1_);
    rx2 = REAL (x2_);
    ry2 = REAL (y2_);
    rout = REAL (out);

    for (size_t i = 0; i < n1; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        double siny1 = sin (ry1 [i] * M_PI / 180.0); // y-value of x data
        double cosy1 = cos (ry1 [i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < n2; j++)
        {
            double siny2 = sin (ry2 [j] * M_PI / 180.0);
            double cosy2 = cos (ry2 [j] * M_PI / 180.0);
            rout [i * n2 + j] = one_vincenty (rx1 [i], ry1 [i],
                    rx2 [j], ry2 [j], siny1, cosy1, siny2, cosy2);
        }
    }

    UNPROTECT (1);

    return out;
}

//' R_cheap_xy_vec
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_xy_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_)
{
    size_t n1 = length (x1_);
    size_t n2 = length (x2_);
    size_t nm = n1 * n2;
    SEXP out = PROTECT (allocVector (REALSXP, nm));
    double *rx1, *ry1, *rx2, *ry2, *rout;
    rx1 = REAL (x1_);
    ry1 = REAL (y1_);
    rx2 = REAL (x2_);
    ry2 = REAL (y2_);
    rout = REAL (out);

    // Get maximal latitude range
    double ymin = 9999.9, ymax = -9999.9;
    for (size_t i = 0; i < n1; i++)
    {
        if (ry1 [i] < ymin)
            ymin = ry1 [i];
        if (ry1 [i] > ymax)
            ymax = ry1 [i];
    }
    for (size_t i = 0; i < n2; i++)
    {
        if (ry2 [i] < ymin)
            ymin = ry2 [i];
        if (ry2 [i] > ymax)
            ymax = ry2 [i];
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    double cosy = cos ((ymin + ymax) / 2.0);

    for (size_t i = 0; i < n1; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = 0; j < n2; j++)
        {
            rout [i * n2 + j] = one_cheap (rx1 [i], ry1 [i],
                    rx2 [j], ry2 [j], cosy);
        }
    }

    UNPROTECT (1);

    return out;
}


//' R_geodesic_xy_vec
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_xy_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_)
{
    size_t n1 = length (x1_);
    size_t n2 = length (x2_);
    size_t nm = n1 * n2;
    SEXP out = PROTECT (allocVector (REALSXP, nm));
    double *rx1, *ry1, *rx2, *ry2, *rout;
    rx1 = REAL (x1_);
    ry1 = REAL (y1_);
    rx2 = REAL (x2_);
    ry2 = REAL (y2_);
    rout = REAL (out);

    for (size_t i = 0; i < n1; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = 0; j < n2; j++)
        {
            rout [i * n2 + j] = one_geodesic (rx1 [i], ry1 [i],
                    rx2 [j], ry2 [j]);
        }
    }

    UNPROTECT (1);

    return out;
}
