#include <math.h>

#include "common.h"

// All values in metres; see https://github.com/hypertidy/geodist/issues/7
// meridian is only used in mapbox cheap distances
const double earth = 6378137.0; // WSG-84 definition
const double meridian = 20003930.0; // length of prime meridian
const double equator = 40007862.917; // IUGG standard equatorial circumference

// All values in metres; see https://github.com/hypertidy/geodist/issues/7
//const double earth = 6378137.0; // WSG-84 definition
//const double meridian = 20003930.0;; // length of prime meridian
//const double equator = 40007862.917; // IUGG standard equatorial circumference

// Core calculations for a single distance measure

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
