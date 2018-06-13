#include <R.h>
#include <Rinternals.h>

#include "dists.h"

//' Haversine for variable x and y
//'
//' @return single distance
//'
//' @note The sxd and syd values could be calculated in arrays, each value of
//' which could be determined with only n operations, rather than the n2 used
//' here. Doing so, however, requires very large C arrays which are often
//' problematic, so this is safer.
//'
//' @noRd
double one_haversine (double x1, double y1, double x2, double y2,
        double cosy1, double cosy2)
{
    double sxd = sin ((x2 - x1) * M_PI / 360.0);
    double syd = sin ((y2 - y1) * M_PI / 360.0);
    double d = syd * syd + cosy1 * cosy2 * sxd * sxd;
    d = 2.0 * earth * asin (sqrt (d));
    return (d);
}

//' Vincenty great circle distance
//' @return single distance
//' @noRd
double one_vincenty (double x1, double y1, double x2, double y2,
        double siny1, double cosy1, double siny2, double cosy2)
{
    double xd = (x2 - x1) * M_PI / 180.0;
    double cxd = cos (xd);
    double n1 = cosy2 * sin (xd); // first term of numerator
    double n2 = cosy1 * siny2 - siny1 * cosy2 * cxd;
    double numerator = n1 * n1 + n2 * n2;
    double denominator = siny1 * siny2 + cosy1 * cosy2 * cxd;
    double d = earth * atan2 (sqrt (numerator), denominator);
    return d;
}
//' mapbox cheap ruler
//' https://blog.mapbox.com/fast-geodesic-approximations-with-cheap-ruler-106f229ad016
double one_cheap (double x1, double y1, double x2, double y2, double cosy)
{
    double dy = meridian * (y1 - y2) / 180.0;
    double dx = equator * (x1 - x2) * cosy / 360.0;
    double d = sqrt (dx * dx + dy * dy);
    return d;
}

//' R_haversine
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    size_t n2 = n * n;
    //Rprintf ("n = %d ; len = %d \n", n, n2);
    SEXP out = PROTECT (allocVector (REALSXP, n2));
    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    double cosy1 [n]; // y-values are indexed in [n+1:n]
    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (rx [n + i] * M_PI / 180.0);
        rout [i * n + i] = 0.0;
    }

    for (size_t i = 0; i < (n - 1); i++)
        for (size_t j = (i + 1); j < n; j++)
        {
            size_t indx1 = i * n + j;
            size_t indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_haversine (rx [i], rx [n + i],
                    rx [j], rx [n + j], cosy1 [i], cosy1 [j]);
        }

    UNPROTECT (1);

    return out;
}

//' R_vincenty
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    size_t n2 = n * n;
    //Rprintf ("n = %d ; len = %d \n", n, n2);
    SEXP out = PROTECT (allocVector (REALSXP, n2));
    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    double siny1 [n], cosy1 [n]; // y-values are indexed in [n+1:n]
    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (rx [n + i] * M_PI / 180.0);
        siny1 [i] = sin (rx [n + i] * M_PI / 180.0);
        rout [i * n + i] = 0.0;
    }

    for (size_t i = 0; i < (n - 1); i++)
        for (size_t j = (i + 1); j < n; j++)
        {
            size_t indx1 = i * n + j;
            size_t indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_vincenty (rx [i], rx [n + i],
                    rx [j], rx [n + j],
                    siny1 [i], cosy1 [i], siny1 [j], cosy1 [j]);
        }

    UNPROTECT (1);

    return out;
}

//' R_cheap
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap (SEXP x_)
{
    size_t n = floor (length (x_) / 2);
    size_t n2 = n * n;
    SEXP out = PROTECT (allocVector (REALSXP, n2));
    double *rx, *rout;
    rx = REAL (x_);
    rout = REAL (out);

    // Get maximal latitude range
    double ymin = 9999.9, ymax = -9999.9;
    for (size_t i = 0; i < n; i++)
    {
        if (rx [n + i] < ymin)
            ymin = rx [n + i];
        else if (rx [n + i] > ymax)
            ymax = rx [n + i];
        rout [i * n + i] = 0.0;
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    double cosy = cos (ymin + ymax) / 2.0;

    for (size_t i = 0; i < (n - 1); i++)
        for (size_t j = (i + 1); j < n; j++)
        {
            size_t indx1 = i * n + j;
            size_t indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_cheap (rx [i], rx [n + i],
                    rx [j], rx [n + j], cosy);
        }

    UNPROTECT (1);

    return out;
}

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
    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < nx; i++)
    {
        double cosy1 = cos (rx [nx + i] * M_PI / 180.0); // y-value of x data
        for (size_t j = 0; j < ny; j++)
        {
            double cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            rout [i * ny + j] = one_haversine (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], cosy1, cosy2);
        }
    }

    UNPROTECT (1);

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
    //Rprintf ("(nx, ny) = (%d , %d )\n", nx, ny);
    SEXP out = PROTECT (allocVector (REALSXP, n2));
    double *rx, *ry, *rout;
    rx = REAL (x_);
    ry = REAL (y_);
    rout = REAL (out);

    for (size_t i = 0; i < nx; i++)
    {
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

    UNPROTECT (1);

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
        else if (rx [nx + i] > ymax)
            ymax = rx [nx + i];
    }
    for (size_t i = 0; i < ny; i++)
    {
        if (ry [ny + i] < ymin)
            ymin = ry [ny + i];
        else if (ry [ny + i] > ymax)
            ymax = ry [ny + i];
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    double cosy = cos (ymin + ymax) / 2.0;

    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            rout [i * ny + j] = one_cheap (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], cosy);
        }
    }

    UNPROTECT (1);

    return out;
}

