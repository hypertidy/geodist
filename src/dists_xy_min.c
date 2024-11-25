#include "dists_xy_min.h"

//' R_haversine_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_xy_min (SEXP x_, SEXP y_)
{
    size_t nx = (size_t) (floor (length (x_) / 2));
    size_t ny = (size_t) (floor (length (y_) / 2));
    
    double *rx, *ry;
    double cosy1, cosy2;
    double d, dmin = 0;
    int *iout;
    int jmin = -1;

    SEXP out = PROTECT (allocVector (INTSXP, nx));

    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);

    iout = INTEGER (out);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
                                     //
        cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data
                                                  //
        dmin = DBL_MAX;
        for (size_t j = 0; j < ny; j++)
        {
            cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            d = one_haversine (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], cosy1, cosy2);
            if (d < dmin)
            {
                dmin = d;
                jmin = j;
            }
        }
        iout [i] = jmin + 1L;
    }

    UNPROTECT (3);

    return out;
}

//' R_vincenty_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty_xy_min (SEXP x_, SEXP y_)
{
    size_t nx = (size_t) (floor (length (x_) / 2));
    size_t ny = (size_t) (floor (length (y_) / 2));

    double *rx, *ry;
    double siny1, cosy1, siny2, cosy2;
    double d, dmin = 0;
    int *iout;
    int jmin = -1;

    SEXP out = PROTECT (allocVector (INTSXP, nx));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    iout = INTEGER (out);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        siny1 = sin (rx [nx + i] * M_PI / 180.0); // y-value of x data
        cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data

        dmin = DBL_MAX;
        for (size_t j = 0; j < ny; j++)
        {
            siny2 = sin (ry [ny + j] * M_PI / 180.0);
            cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            d = one_vincenty (rx [i], ry [j],
                    siny1, cosy1, siny2, cosy2);
            if (d < dmin)
            {
                dmin = d;
                jmin = j;
            }
        }
        iout [i] = jmin + 1L;
    }

    UNPROTECT (3);

    return out;
}

//' R_cheap_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap_xy_min (SEXP x_, SEXP y_)
{
    size_t nx = (size_t) (floor (length (x_) / 2));
    size_t ny = (size_t) (floor (length (y_) / 2));

    double *rx, *ry;
    double ymin = 9999.9, ymax = -9999.9;
    double cosy;
    double d, dmin = 0;
    int *iout;
    int jmin = -1;

    SEXP out = PROTECT (allocVector (INTSXP, nx));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    iout = INTEGER (out);

    // Get maximal latitude range
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
    cosy = cos ((ymin + ymax) / 2.0);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        dmin = DBL_MAX;
        for (size_t j = 0; j < ny; j++)
        {
            d = one_cheap (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], cosy);
            if (d < dmin)
            {
                dmin = d;
                jmin = j;
            }
        }
        iout [i] = jmin + 1L;
    }

    UNPROTECT (3);

    return out;
}


//' R_geodesic_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic_xy_min (SEXP x_, SEXP y_)
{
    size_t nx = (size_t) (floor (length (x_) / 2));
    size_t ny = (size_t) (floor (length (y_) / 2));

    double *rx, *ry;
    double d, dmin = 0;
    int *iout;
    int jmin = -1;

    SEXP out = PROTECT (allocVector (INTSXP, nx));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    iout = INTEGER (out);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        dmin = DBL_MAX;
        for (size_t j = 0; j < ny; j++)
        {
            d = one_geodesic (rx [i], rx [nx + i],
                    ry [j], ry [ny + j]);
            if (d < dmin)
            {
                dmin = d;
                jmin = j;
            }
        }
        iout [i] = jmin + 1L;
    }

    UNPROTECT (3);

    return out;
}
