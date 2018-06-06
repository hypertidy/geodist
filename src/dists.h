#include <vector>
#include <cmath>

#include <Rcpp.h>

const double earth = 6378.137; // value used in geosphere::distHaversine

constexpr float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
constexpr int INFINITE_INT =  std::numeric_limits<int>::max ();

double one_haversine (double x1, double y1, double x2, double y2,
        double cosy1, double cosy2);

Rcpp::NumericMatrix rcpp_haversine (Rcpp::NumericMatrix x);
Rcpp::NumericMatrix rcpp_haversine_xy (Rcpp::NumericMatrix x,
        Rcpp::NumericMatrix y);
