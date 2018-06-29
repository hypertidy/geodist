#include <math.h>

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

//' Vincenty ellipsoidal distance
//' @return single distance
//' @noRd
double one_vincenty_ellips (double U1, double U2, double L)
{
    const double tol = 1.0e-12;
    struct VincEllips ve;
    ve.L = L;
    double lambda = L;
    double lambda_old = one_lambda (U1, U2, &ve, lambda),
           delta = 1.0;
    int count = 0, check = 0;;
    while (delta > tol)
    {
        double lambda = one_lambda (U1, U2, &ve, lambda_old);
        delta = fabs (lambda - lambda_old);
        lambda_old = lambda;
        count++;
        if (count > 1e3)
        {
            check = 1;
            break;
        }
    }
    if (check > 0)
        Rprintf ("ERROR: Vincenty ellipsoid failed to converge\n");

    double b2 = b * b,
           u2 = ve.cos2_alpha * ((earth * earth - b2) / b2),
           A = 1.0 + (u2 / 16384.0) * (4096.0 + u2 * (-768.0 +
                       u2 * (320.0 - 175 * u2))),
           B = (u2 / 1024.0) * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)));
    double delta_sig = B * ve.sin_sigma * (ve.cos_2sig_m + 0.25 * B *
            (ve.cos_sigma * (-1.0 + 2.0 * ve.cos_2sig_m) - (B / 6.0) *
             ve.cos_2sig_m * (-3.0 + 4.0 * ve.sin_sigma * ve.sin_sigma) *
             (-3.0 + 4.0 * ve.cos_2sig_m)));
    double s = b * A * (ve.sigma - delta_sig);

    return s;
}

double one_lambda (double U1, double U2, struct VincEllips *ve, double lambda)
{
    double sU1 = sin (U1), sU2 = sin (U2),
           cU1 = cos (U1), cU2 = cos (U2),
           sL = sin (lambda), cL = cos (lambda);
    double t1 = cU2 * sL,
           t2 = cU1 * sU2,
           t3 = sU1 * cU2 * cL;

    struct VincEllips ve_temp;

    ve_temp.sin_sigma = sqrt (t1 * t1 + (t2 - t3) * (t2 - t3));
    ve_temp.cos_sigma = sU1 * sU2 + cU1 * cU2 * cL;
    ve_temp.sigma = atan2 (ve_temp.sin_sigma, ve_temp.cos_sigma);

    ve_temp.sin_alpha = cU1 * cU2 * sL / ve_temp.sin_sigma;
    ve_temp.cos2_alpha = 1 - ve_temp.sin_alpha * ve_temp.sin_alpha;
    ve_temp.cos_2sig_m = ve_temp.cos_sigma - 2 * sU1 * sU2 / ve_temp.cos2_alpha;

    double C = (flattening / 16.0) * ve_temp.cos2_alpha *
        (4.0 + flattening * (4.0 - 3.0 * ve_temp.cos2_alpha));

    double lambda_new = ve_temp.L + (1.0 - C) * flattening * ve_temp.sin_alpha *
        (ve_temp.sigma + C * ve_temp.sin_sigma * (ve_temp.cos_2sig_m +
                                                  C * ve_temp.cos_sigma *
                                          (-1.0 + 2.0 * ve_temp.cos_2sig_m)));
    *ve = ve_temp;

    return lambda_new;
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
