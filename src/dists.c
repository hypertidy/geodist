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
