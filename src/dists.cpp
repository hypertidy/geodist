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
    double d = sin (yd / 2.0) * sin (yd / 2.0) + cos (y2 * M_PI / 180.0) *
        cos (y1 * M_PI / 180.0) * sin (xd / 2.0) * sin (xd / 2.0);
    d = 2.0 * earth * asin (sqrt (d));
    return (d);
}
