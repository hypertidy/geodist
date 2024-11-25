#include "dists_xy.h"

//' R_haversine_xy
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @param y_ Additional vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine_xy (SEXP x_, SEXP y_)
{
    size_t nx = (size_t) (floor (length (x_) / 2));
    size_t ny = (size_t) (floor (length (y_) / 2));
    size_t n2 = nx * ny;
    //Rprintf ("(nx, ny) = (%d , %d )\n", nx, ny);
    //
    double *rx, *ry, *rout;
    double cosy1, cosy2;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < ny; j++)
        {
            cosy2 = cos (ry [ny + j] * M_PI / 180.0);
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
    size_t nx = (size_t) (floor (length (x_) / 2));
    size_t ny = (size_t) (floor (length (y_) / 2));
    size_t n2 = nx * ny;

    double *rx, *ry, *rout;
    double siny1, cosy1, siny2, cosy2;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < nx; i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        siny1 = sin (rx [nx + i] * M_PI / 180.0); // y-value of x data
        cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < ny; j++)
        {
            siny2 = sin (ry [ny + j] * M_PI / 180.0);
            cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            rout [i * ny + j] = one_vincenty (rx [i], ry [j],
                    siny1, cosy1, siny2, cosy2);
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
    size_t nx = (size_t) (floor (length (x_) / 2));
    size_t ny = (size_t) (floor (length (y_) / 2));
    size_t n2 = nx * ny;

    double *rx, *ry, *rout;
    double ymin = 9999.9, ymax = -9999.9;
    double cosy;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

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
    size_t nx = (size_t) (floor (length (x_) / 2));
    size_t ny = (size_t) (floor (length (y_) / 2));
    size_t n2 = nx * ny;

    double *rx, *ry, *rout;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));
    y_ = PROTECT (Rf_coerceVector (y_, REALSXP));

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
