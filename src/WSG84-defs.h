#ifndef WSG_H
#define WSG_H

// All values in metres; see https://github.com/hypertidy/geodist/issues/7
// meridian is only used in mapbox cheap distances
static const double earth = 6378137.0; // WSG-84 definition
static const double meridian = 20003930.0; // length of prime meridian
static const double equator = 40007862.917; // IUGG standard equatorial circumference
static const double flattening = 1.0 / 298.257223563; // flattening of ellipsoid
//static const double b = (1.0 - flattening) * earth; // only for vincenty ellips
static const double b = 6356752.314245;

#endif /* COMMON_H */
