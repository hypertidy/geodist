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

    std::vector <double> cosy1 (static_cast <size_t> (nx));
    for (size_t i = 0; i < nx; i++)
        cosy1 [i] = cos (x (i, 1) * M_PI / 180.0);

    for (size_t i = 0; i < (nx - 1); i++)
        for (size_t j = (i + 1); j < nx; j++)
        {
            res (i, j) = res (j, i) =
                one_haversine (x (i, 0), x (i, 1), x (j, 0), x (j, 1),
                        cosy1 [i], cosy1 [j]);
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
    {
        double cosy1 = cos (x (i, 1) * M_PI / 180.0);

        for (size_t j = 0; j < ny; j++)
        {
            double cosy2 = cos (y (j, 1) * M_PI / 180.0);
            res (i, j) = one_haversine (x (i, 0), x (i, 1), y (j, 0), y (j, 1),
                cosy1, cosy2);
        }
    }
    return res;
}
