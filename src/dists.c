#include <R.h>
#include <Rinternals.h>

#include "dists.h"

//' haversine
//'
//' Haversine for variable x and y
//'
//' @return single distance
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
            double siny2 = cos (ry [ny + j] * M_PI / 180.0);
            double cosy2 = cos (ry [ny + j] * M_PI / 180.0);
            rout [i * ny + j] = one_vincenty (rx [i], rx [nx + i],
                    ry [j], ry [ny + j], siny1, cosy1, siny2, cosy2);
        }
    }

    UNPROTECT (1);

    return out;
}
