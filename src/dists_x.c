#include "dists_x.h"

//' R_haversine
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_haversine (SEXP x_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    //Rprintf ("n = %d ; len = %d \n", n, n2);
    size_t n2 = n * n;
    double cosy1 [n]; // y-values are indexed in [n+1:n]

    double *rx, *rout;
    size_t indx1, indx2;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    rx = REAL (x_);
    rout = REAL (out);

    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (rx [n + i] * M_PI / 180.0);
        rout [i * n + i] = 0.0;
    }

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = (i + 1); j < n; j++)
        {
            indx1 = i * n + j;
            indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_haversine (rx [i], rx [n + i],
                    rx [j], rx [n + j], cosy1 [i], cosy1 [j]);
        }
    }

    UNPROTECT (2);

    return out;
}

//' R_vincenty
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_vincenty (SEXP x_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    size_t n2 = n * n;
    //Rprintf ("n = %d ; len = %d \n", n, n2);

    double *rx, *rout;
    double siny1 [n], cosy1 [n]; // y-values are indexed in [n+1:n]
    size_t indx1, indx2;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    rx = REAL (x_);
    rout = REAL (out);

    for (size_t i = 0; i < n; i++)
    {
        cosy1 [i] = cos (rx [n + i] * M_PI / 180.0);
        siny1 [i] = sin (rx [n + i] * M_PI / 180.0);
        rout [i * n + i] = 0.0;
    }

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = (i + 1); j < n; j++)
        {
            indx1 = i * n + j;
            indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_vincenty (rx [i], rx [j],
                    siny1 [i], cosy1 [i], siny1 [j], cosy1 [j]);
        }
    }

    UNPROTECT (2);

    return out;
}

//' R_cheap
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_cheap (SEXP x_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    size_t n2 = n * n;

    double *rx, *rout;
    double ymin = 9999.9, ymax = -9999.9;
    double cosy;
    size_t indx1, indx2;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    rx = REAL (x_);
    rout = REAL (out);

    // Get maximal latitude range
    for (size_t i = 0; i < n; i++)
    {
        if (rx [n + i] < ymin)
            ymin = rx [n + i];
        if (rx [n + i] > ymax)
            ymax = rx [n + i];
        rout [i * n + i] = 0.0;
    }
    // and set constant cosine multiplier
    ymin = ymin * M_PI / 180;
    ymax = ymax * M_PI / 180;
    cosy = cos ((ymin + ymax) / 2.0);

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = (i + 1); j < n; j++)
        {
            indx1 = i * n + j;
            indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_cheap (rx [i], rx [n + i],
                    rx [j], rx [n + j], cosy);
        }
    }

    UNPROTECT (2);

    return out;
}

//' R_geodesic
//' @param x_ Single vector of x-values in [1:n], y-values in [n+(1:n)]
//' @noRd
SEXP R_geodesic (SEXP x_)
{
    size_t n = (size_t) (floor (length (x_) / 2));
    size_t n2 = n * n;
    double *rx, *rout;
    size_t indx1, indx2;

    SEXP out = PROTECT (allocVector (REALSXP, n2));
    x_ = PROTECT (Rf_coerceVector (x_, REALSXP));

    rx = REAL (x_);
    rout = REAL (out);

    for (size_t i = 0; i < n; i++)
        rout [i * n + i] = 0.0;

    for (size_t i = 0; i < (n - 1); i++)
    {
        if (i % 1000 == 0)
            R_CheckUserInterrupt (); // # nocov
        for (size_t j = (i + 1); j < n; j++)
        {
            indx1 = i * n + j;
            indx2 = j * n + i;
            rout [indx1] = rout [indx2] = one_geodesic (rx [i], rx [n + i],
                    rx [j], rx [n + j]);
        }
    }

    UNPROTECT (2);

    return out;
}
