#include "dists.h"

//' haversine
//'
//' Haversine for variable x and y
//'
//' @return single distance
//'
//' @noRd
double one_haversine (double x1, double y1, double x2, double y2)
{
    double xd = (x2 - x1) * M_PI / 180.0;
    double yd = (y2 - y1) * M_PI / 180.0;
    double sxd = sin (xd / 2.0);
    double syd = sin (yd / 2.0);
    double d = syd * syd + cos (y2 * M_PI / 180.0) *
        cos (y1 * M_PI / 180.0) * sxd * sxd;
    d = 2.0 * earth * asin (sqrt (d));
    return (d);
}

//' rcpp_haversine
//'
//' @return single distance
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_haversine (Rcpp::NumericMatrix x)
{
    const int nx = x.nrow ();
    Rcpp::NumericMatrix res (nx, nx);
    for (size_t i = 0; i < nx; i++)
        res (i, i) = 0.0;
    for (size_t i = 0; i < (nx - 1); i++)
        for (size_t j = (i + 1); j < nx; j++)
        {
            res (i, j) = res (j, i) =
                one_haversine (x (i, 0), x (i, 1), x (j, 0), x (j, 1));
        }
    return res;
}

//' rcpp_haversine
//'
//' @return single distance
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_haversine_xy (Rcpp::NumericMatrix x,
        Rcpp::NumericMatrix y)
{
    const int nx = x.nrow (),
              ny = y.nrow ();
    Rcpp::NumericMatrix res (nx, ny);
    for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
        {
            res (i, j) = one_haversine (x (i, 0), x (i, 1), x (j, 0), x (j, 1));
        }
    return res;
}
