#ifndef COMMON_H
#define COMMON_H

double one_haversine (double x1, double y1, double x2, double y2,
        double cosy1, double cosy2);
double one_vincenty (double x1, double y1, double x2, double y2,
        double siny1, double cosy1, double siny2, double cosy2);

double one_cheap (double x1, double y1, double x2, double y2, double cosy);
double one_geodesic (double x1, double y1, double x2, double y2);

#endif /* COMMON_H */
