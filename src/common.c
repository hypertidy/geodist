#include <math.h>
#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"
#include "geodesic.h"


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
double one_vincenty (double x1, double x2,
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

//' Karney (2013) geodesic
//' https://geographiclib.sourceforge.io/geod.html
//' https://link.springer.com/content/pdf/10.1007/s00190-012-0578-z.pdf
double one_geodesic (double x1, double y1, double x2, double y2)
{
    //double lat1, lon1, azi1, lat2, lon2, azi2, s12;
    double azi1, azi2, s12;
    struct geod_geodesic g;

    geod_init(&g, earth, flattening);
    geod_inverse(&g, y1, x1, y2, x2, &s12, &azi1, &azi2);
    return s12;
}
