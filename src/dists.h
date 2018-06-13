const double earth = 6378.137; // value used in geosphere::distHaversine
const double meridian = 12429.9; // length of prime meridian in metres
const double equator = 40007.862917; // equatorial circumference

double one_haversine (double x1, double y1, double x2, double y2,
        double cosy1, double cosy2);
double one_vincenty (double x1, double y1, double x2, double y2,
        double siny1, double cosy1, double siny2, double cosy2);
double one_cheap (double x1, double y1, double x2, double y2, double cosy);
