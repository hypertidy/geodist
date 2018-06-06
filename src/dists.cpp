#include "dists.h"

//' rcpp_haversine
//'
//' @return single distance
//'
//' @noRd
// [[Rcpp::export]]
double rcpp_haversine (double x1, double y1, double x2, double y2)
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
