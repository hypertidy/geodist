/**
 * \file geodesic.c
 * \brief Implementation of the geodesic routines in C
 *
 * For the full documentation see geodesic.h.
 **********************************************************************/

/** @cond SKIP */

/*
 * This is a C implementation of the geodesic algorithms described in
 *
 *   C. F. F. Karney,
 *   Algorithms for geodesics,
 *   J. Geodesy <b>87</b>, 43--55 (2013);
 *   https://doi.org/10.1007/s00190-012-0578-z
 *   Addenda: https://geographiclib.sourceforge.io/geod-addenda.html
 *
 * See the comments in geodesic.h for documentation.
 *
 * Copyright (c) Charles Karney (2012-2017) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 */

#include "geodesic.h"
#include <math.h>

#if !defined(HAVE_C99_MATH)
#define HAVE_C99_MATH 0
#endif

#define GEOGRAPHICLIB_GEODESIC_ORDER 6
#define nA1   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC1   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC1p  GEOGRAPHICLIB_GEODESIC_ORDER
#define nA2   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC2   GEOGRAPHICLIB_GEODESIC_ORDER
#define nA3   GEOGRAPHICLIB_GEODESIC_ORDER
#define nA3x  nA3
#define nC3   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC3x  ((nC3 * (nC3 - 1)) / 2)
#define nC4   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC4x  ((nC4 * (nC4 + 1)) / 2)
#define nC    (GEOGRAPHICLIB_GEODESIC_ORDER + 1)

typedef double real;
typedef int boolx;

static unsigned init = 0;
static const int FALSE = 0;
static const int TRUE = 1;
static unsigned digits, maxit1, maxit2;
static real epsilon, realmin, pi, degree, NaN,
  tiny, tol0, tol1, tol2, tolb, xthresh;

static void Init() {
  if (!init) {
#if defined(__DBL_MANT_DIG__)
    digits = __DBL_MANT_DIG__;
#else
    digits = 53;
#endif
#if defined(__DBL_EPSILON__)
    epsilon = __DBL_EPSILON__;
#else
    epsilon = pow(0.5, digits - 1);
#endif
#if defined(__DBL_MIN__)
    realmin = __DBL_MIN__;
#else
    realmin = pow(0.5, 1022);
#endif
#if defined(M_PI)
    pi = M_PI;
#else
    pi = atan2(0.0, -1.0);
#endif
    maxit1 = 20;
    maxit2 = maxit1 + digits + 10;
    tiny = sqrt(realmin);
    tol0 = epsilon;
    /* Increase multiplier in defn of tol1 from 100 to 200 to fix inverse case
     * 52.784459512564 0 -52.784459512563990912 179.634407464943777557
     * which otherwise failed for Visual Studio 10 (Release and Debug) */
    tol1 = 200 * tol0;
    tol2 = sqrt(tol0);
    /* Check on bisection interval */
    tolb = tol0 * tol2;
    xthresh = 1000 * tol2;
    degree = pi/180;
    {
      real minus1 = -1;
      NaN = sqrt(minus1);
    }
    init = 1;
  }
}

enum captype {
  CAP_NONE = 0U,
  CAP_C1   = 1U<<0,
  CAP_C1p  = 1U<<1,
  CAP_C2   = 1U<<2,
  CAP_C3   = 1U<<3,
  CAP_C4   = 1U<<4,
  CAP_ALL  = 0x1FU,
  OUT_ALL  = 0x7F80U
};

static real sq(real x) { return x * x; }
#if HAVE_C99_MATH
#define atanhx atanh
#define copysignx copysign
#define hypotx hypot
#define cbrtx cbrt
#else
static real log1px(real x) {
  volatile real
    y = 1 + x,
    z = y - 1;
  /* Here's the explanation for this magic: y = 1 + z, exactly, and z
   * approx x, thus log(y)/z (which is nearly constant near z = 0) returns
   * a good approximation to the true log(1 + x)/x.  The multiplication x *
   * (log(y)/z) introduces little additional error. */
  return z == 0 ? x : x * log(y) / z;
}

static real atanhx(real x) {
  real y = fabs(x);             /* Enforce odd parity */
  y = log1px(2 * y/(1 - y))/2;
  return x < 0 ? -y : y;
}

static real copysignx(real x, real y) {
  return fabs(x) * (y < 0 || (y == 0 && 1/y < 0) ? -1 : 1);
}

static real hypotx(real x, real y)
{ return sqrt(x * x + y * y); }

static real cbrtx(real x) {
  real y = pow(fabs(x), 1/(real)(3)); /* Return the real cube root */
  return x < 0 ? -y : y;
}
#endif

static real sumx(real u, real v, real* t) {
  volatile real s = u + v;
  volatile real up = s - v;
  volatile real vpp = s - up;
  up -= u;
  vpp -= v;
  if (t) *t = -(up + vpp);
  /* error-free sum:
   * u + v =       s      + t
   *       = round(u + v) + t */
  return s;
}

static real polyval(int N, const real p[], real x) {
  real y = N < 0 ? 0 : *p++;
  while (--N >= 0) y = y * x + *p++;
  return y;
}

/* mimic C++ std::min and std::max */
static real minx(real a, real b)
{ return (b < a) ? b : a; }

static real maxx(real a, real b)
{ return (a < b) ? b : a; }

static void swapx(real* x, real* y)
{ real t = *x; *x = *y; *y = t; }

static void norm2(real* sinx, real* cosx) {
  real r = hypotx(*sinx, *cosx);
  *sinx /= r;
  *cosx /= r;
}

static real AngNormalize(real x) {
#if HAVE_C99_MATH
  x = remainder(x, (real)(360));
  return x != -180 ? x : 180;
#else
  x = fmod(x, (real)(360));
  return x <= -180 ? x + 360 : (x <= 180 ? x : x - 360);
#endif
}

static real LatFix(real x)
{ return fabs(x) > 90 ? NaN : x; }

static real AngDiff(real x, real y, real* e) {
  real t, d = AngNormalize(sumx(AngNormalize(-x), AngNormalize(y), &t));
  /* Here y - x = d + t (mod 360), exactly, where d is in (-180,180] and
   * abs(t) <= eps (eps = 2^-45 for doubles).  The only case where the
   * addition of t takes the result outside the range (-180,180] is d = 180
   * and t > 0.  The case, d = -180 + eps, t = -eps, can't happen, since
   * sum would have returned the exact result in such a case (i.e., given t
   * = 0). */
  return sumx(d == 180 && t > 0 ? -180 : d, t, e);
}

static real AngRound(real x) {
  const real z = 1/(real)(16);
  volatile real y;
  if (x == 0) return 0;
  y = fabs(x);
  /* The compiler mustn't "simplify" z - (z - y) to y */
  y = y < z ? z - (z - y) : y;
  return x < 0 ? -y : y;
}

static void sincosdx(real x, real* sinx, real* cosx) {
  /* In order to minimize round-off errors, this function exactly reduces
   * the argument to the range [-45, 45] before converting it to radians. */
  real r, s, c; int q;
#if HAVE_C99_MATH && !defined(__GNUC__)
  /* Disable for gcc because of bug in glibc version < 2.22, see
   * https://sourceware.org/bugzilla/show_bug.cgi?id=17569 */
  r = remquo(x, (real)(90), &q);
#else
  r = fmod(x, (real)(360));
  q = (int)(floor(r / 90 + (real)(0.5)));
  r -= 90 * q;
#endif
  /* now abs(r) <= 45 */
  r *= degree;
  /* Possibly could call the gnu extension sincos */
  s = sin(r); c = cos(r);
  switch ((unsigned)q & 3U) {
  case 0U: *sinx =  s; *cosx =  c; break;
  case 1U: *sinx =  c; *cosx = -s; break;
  case 2U: *sinx = -s; *cosx = -c; break;
  default: *sinx = -c; *cosx =  s; break; /* case 3U */
  }
  if (x != 0) { *sinx += (real)(0); *cosx += (real)(0); }
}

static real atan2dx(real y, real x) {
  /* In order to minimize round-off errors, this function rearranges the
   * arguments so that result of atan2 is in the range [-pi/4, pi/4] before
   * converting it to degrees and mapping the result to the correct
   * quadrant. */
  int q = 0; real ang;
  if (fabs(y) > fabs(x)) { swapx(&x, &y); q = 2; }
  if (x < 0) { x = -x; ++q; }
  /* here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4] */
  ang = atan2(y, x) / degree;
  switch (q) {
    /* Note that atan2d(-0.0, 1.0) will return -0.  However, we expect that
     * atan2d will not be called with y = -0.  If need be, include
     *
     *   case 0: ang = 0 + ang; break;
     */
  case 1: ang = (y >= 0 ? 180 : -180) - ang; break;
  case 2: ang =  90 - ang; break;
  case 3: ang = -90 + ang; break;
  }
  return ang;
}

static void A3coeff(struct geod_geodesic* g);
static void C3coeff(struct geod_geodesic* g);
static void C4coeff(struct geod_geodesic* g);
static real SinCosSeries(boolx sinp,
                         real sinx, real cosx,
                         const real c[], int n);
static void Lengths(const struct geod_geodesic* g,
                    real eps, real sig12,
                    real ssig1, real csig1, real dn1,
                    real ssig2, real csig2, real dn2,
                    real cbet1, real cbet2,
                    real* ps12b, real* pm12b, real* pm0,
                    real* pM12, real* pM21,
                    /* Scratch area of the right size */
                    real Ca[]);
static real Astroid(real x, real y);
static real InverseStart(const struct geod_geodesic* g,
                         real sbet1, real cbet1, real dn1,
                         real sbet2, real cbet2, real dn2,
                         real lam12, real slam12, real clam12,
                         real* psalp1, real* pcalp1,
                         /* Only updated if return val >= 0 */
                         real* psalp2, real* pcalp2,
                         /* Only updated for short lines */
                         real* pdnm,
                         /* Scratch area of the right size */
                         real Ca[]);
static real Lambda12(const struct geod_geodesic* g,
                     real sbet1, real cbet1, real dn1,
                     real sbet2, real cbet2, real dn2,
                     real salp1, real calp1,
                     real slam120, real clam120,
                     real* psalp2, real* pcalp2,
                     real* psig12,
                     real* pssig1, real* pcsig1,
                     real* pssig2, real* pcsig2,
                     real* peps,
                     real* pgomg12,
                     boolx diffp, real* pdlam12,
                     /* Scratch area of the right size */
                     real Ca[]);
static real A3f(const struct geod_geodesic* g, real eps);
static void C3f(const struct geod_geodesic* g, real eps, real c[]);
static void C4f(const struct geod_geodesic* g, real eps, real c[]);
static real A1m1f(real eps);
static void C1f(real eps, real c[]);
static void C1pf(real eps, real c[]);
static real A2m1f(real eps);
static void C2f(real eps, real c[]);
static int transit(real lon1, real lon2);
static int transitdirect(real lon1, real lon2);
static void accini(real s[]);
static void acccopy(const real s[], real t[]);
static void accadd(real s[], real y);
static real accsum(const real s[], real y);
static void accneg(real s[]);

void geod_init(struct geod_geodesic* g, real a, real f) {
  if (!init) Init();
  g->a = a;
  g->f = f;
  g->f1 = 1 - g->f;
  g->e2 = g->f * (2 - g->f);
  g->ep2 = g->e2 / sq(g->f1);   /* e2 / (1 - e2) */
  g->n = g->f / ( 2 - g->f);
  g->b = g->a * g->f1;
  g->c2 = (sq(g->a) + sq(g->b) *
           (g->e2 == 0 ? 1 :
            (g->e2 > 0 ? atanhx(sqrt(g->e2)) : atan(sqrt(-g->e2))) /
            sqrt(fabs(g->e2))))/2; /* authalic radius squared */
  /* The sig12 threshold for "really short".  Using the auxiliary sphere
   * solution with dnm computed at (bet1 + bet2) / 2, the relative error in the
   * azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.  (Error
   * measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a given f and
   * sig12, the max error occurs for lines near the pole.  If the old rule for
   * computing dnm = (dn1 + dn2)/2 is used, then the error increases by a
   * factor of 2.)  Setting this equal to epsilon gives sig12 = etol2.  Here
   * 0.1 is a safety factor (error decreased by 100) and max(0.001, abs(f))
   * stops etol2 getting too large in the nearly spherical case. */
  g->etol2 = 0.1 * tol2 /
    sqrt( maxx((real)(0.001), fabs(g->f)) * minx((real)(1), 1 - g->f/2) / 2 );

  A3coeff(g);
  C3coeff(g);
  C4coeff(g);
}

static real geod_geninverse_int(const struct geod_geodesic* g,
                                real lat1, real lon1, real lat2, real lon2,
                                real* ps12,
                                real* psalp1, real* pcalp1,
                                real* psalp2, real* pcalp2,
                                real* pm12, real* pM12, real* pM21,
                                real* pS12) {
  real s12 = 0, m12 = 0, M12 = 0, M21 = 0, S12 = 0;
  real lon12, lon12s;
  int latsign, lonsign, swapp;
  real sbet1, cbet1, sbet2, cbet2, s12x = 0, m12x = 0;
  real dn1, dn2, lam12, slam12, clam12;
  real a12 = 0, sig12, calp1 = 0, salp1 = 0, calp2 = 0, salp2 = 0;
  real Ca[nC];
  boolx meridian;
  /* somg12 > 1 marks that it needs to be calculated */
  real omg12 = 0, somg12 = 2, comg12 = 0;

  unsigned outmask =
    (ps12 ? GEOD_DISTANCE : 0U) |
    (pm12 ? GEOD_REDUCEDLENGTH : 0U) |
    (pM12 || pM21 ? GEOD_GEODESICSCALE : 0U) |
    (pS12 ? GEOD_AREA : 0U);

  outmask &= OUT_ALL;
  /* Compute longitude difference (AngDiff does this carefully).  Result is
   * in [-180, 180] but -180 is only for west-going geodesics.  180 is for
   * east-going and meridional geodesics. */
  lon12 = AngDiff(lon1, lon2, &lon12s);
  /* Make longitude difference positive. */
  lonsign = lon12 >= 0 ? 1 : -1;
  /* If very close to being on the same half-meridian, then make it so. */
  lon12 = lonsign * AngRound(lon12);
  lon12s = AngRound((180 - lon12) - lonsign * lon12s);
  lam12 = lon12 * degree;
  if (lon12 > 90) {
    sincosdx(lon12s, &slam12, &clam12);
    clam12 = -clam12;
  } else
    sincosdx(lon12, &slam12, &clam12);

  /* If really close to the equator, treat as on equator. */
  lat1 = AngRound(LatFix(lat1));
  lat2 = AngRound(LatFix(lat2));
  /* Swap points so that point with higher (abs) latitude is point 1
   * If one latitude is a nan, then it becomes lat1. */
  swapp = fabs(lat1) < fabs(lat2) ? -1 : 1;
  if (swapp < 0) {
    lonsign *= -1;
    swapx(&lat1, &lat2);
  }
  /* Make lat1 <= 0 */
  latsign = lat1 < 0 ? 1 : -1;
  lat1 *= latsign;
  lat2 *= latsign;
  /* Now we have
   *
   *     0 <= lon12 <= 180
   *     -90 <= lat1 <= 0
   *     lat1 <= lat2 <= -lat1
   *
   * longsign, swapp, latsign register the transformation to bring the
   * coordinates to this canonical form.  In all cases, 1 means no change was
   * made.  We make these transformations so that there are few cases to
   * check, e.g., on verifying quadrants in atan2.  In addition, this
   * enforces some symmetries in the results returned. */

  sincosdx(lat1, &sbet1, &cbet1); sbet1 *= g->f1;
  /* Ensure cbet1 = +epsilon at poles */
  norm2(&sbet1, &cbet1); cbet1 = maxx(tiny, cbet1);

  sincosdx(lat2, &sbet2, &cbet2); sbet2 *= g->f1;
  /* Ensure cbet2 = +epsilon at poles */
  norm2(&sbet2, &cbet2); cbet2 = maxx(tiny, cbet2);

  /* If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
   * |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
   * a better measure.  This logic is used in assigning calp2 in Lambda12.
   * Sometimes these quantities vanish and in that case we force bet2 = +/-
   * bet1 exactly.  An example where is is necessary is the inverse problem
   * 48.522876735459 0 -48.52287673545898293 179.599720456223079643
   * which failed with Visual Studio 10 (Release and Debug) */

  if (cbet1 < -sbet1) {
    if (cbet2 == cbet1)
      sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
  } else {
    if (fabs(sbet2) == -sbet1)
      cbet2 = cbet1;
  }

  dn1 = sqrt(1 + g->ep2 * sq(sbet1));
  dn2 = sqrt(1 + g->ep2 * sq(sbet2));

  meridian = lat1 == -90 || slam12 == 0;

  if (meridian) {

    /* Endpoints are on a single full meridian, so the geodesic might lie on
     * a meridian. */

    real ssig1, csig1, ssig2, csig2;
    calp1 = clam12; salp1 = slam12; /* Head to the target longitude */
    calp2 = 1; salp2 = 0;           /* At the target we're heading north */

    /* tan(bet) = tan(sig) * cos(alp) */
    ssig1 = sbet1; csig1 = calp1 * cbet1;
    ssig2 = sbet2; csig2 = calp2 * cbet2;

    /* sig12 = sig2 - sig1 */
    sig12 = atan2(maxx((real)(0), csig1 * ssig2 - ssig1 * csig2),
                                  csig1 * csig2 + ssig1 * ssig2);
    Lengths(g, g->n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
            cbet1, cbet2, &s12x, &m12x, 0,
            outmask & GEOD_GEODESICSCALE ? &M12 : 0,
            outmask & GEOD_GEODESICSCALE ? &M21 : 0,
            Ca);
    /* Add the check for sig12 since zero length geodesics might yield m12 <
     * 0.  Test case was
     *
     *    echo 20.001 0 20.001 0 | GeodSolve -i
     *
     * In fact, we will have sig12 > pi/2 for meridional geodesic which is
     * not a shortest path. */
    if (sig12 < 1 || m12x >= 0) {
      /* Need at least 2, to handle 90 0 90 180 */
      if (sig12 < 3 * tiny)
        sig12 = m12x = s12x = 0;
      m12x *= g->b;
      s12x *= g->b;
      a12 = sig12 / degree;
    } else
      /* m12 < 0, i.e., prolate and too close to anti-podal */
      meridian = FALSE;
  }

  if (!meridian &&
      sbet1 == 0 &&           /* and sbet2 == 0 */
      /* Mimic the way Lambda12 works with calp1 = 0 */
      (g->f <= 0 || lon12s >= g->f * 180)) {

    /* Geodesic runs along equator */
    calp1 = calp2 = 0; salp1 = salp2 = 1;
    s12x = g->a * lam12;
    sig12 = omg12 = lam12 / g->f1;
    m12x = g->b * sin(sig12);
    if (outmask & GEOD_GEODESICSCALE)
      M12 = M21 = cos(sig12);
    a12 = lon12 / g->f1;

  } else if (!meridian) {

    /* Now point1 and point2 belong within a hemisphere bounded by a
     * meridian and geodesic is neither meridional or equatorial. */

    /* Figure a starting point for Newton's method */
    real dnm = 0;
    sig12 = InverseStart(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                         lam12, slam12, clam12,
                         &salp1, &calp1, &salp2, &calp2, &dnm,
                         Ca);

    if (sig12 >= 0) {
      /* Short lines (InverseStart sets salp2, calp2, dnm) */
      s12x = sig12 * g->b * dnm;
      m12x = sq(dnm) * g->b * sin(sig12 / dnm);
      if (outmask & GEOD_GEODESICSCALE)
        M12 = M21 = cos(sig12 / dnm);
      a12 = sig12 / degree;
      omg12 = lam12 / (g->f1 * dnm);
    } else {

      /* Newton's method.  This is a straightforward solution of f(alp1) =
       * lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
       * root in the interval (0, pi) and its derivative is positive at the
       * root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
       * alp1.  During the course of the iteration, a range (alp1a, alp1b) is
       * maintained which brackets the root and with each evaluation of
       * f(alp) the range is shrunk, if possible.  Newton's method is
       * restarted whenever the derivative of f is negative (because the new
       * value of alp1 is then further from the solution) or if the new
       * estimate of alp1 lies outside (0,pi); in this case, the new starting
       * guess is taken to be (alp1a + alp1b) / 2. */
      real ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0, domg12 = 0;
      unsigned numit = 0;
      /* Bracketing range */
      real salp1a = tiny, calp1a = 1, salp1b = tiny, calp1b = -1;
      boolx tripn, tripb;
      for (tripn = FALSE, tripb = FALSE; numit < maxit2; ++numit) {
        /* the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
         * WGS84 and random input: mean = 2.85, sd = 0.60 */
        real dv = 0,
          v = Lambda12(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                        slam12, clam12,
                        &salp2, &calp2, &sig12, &ssig1, &csig1, &ssig2, &csig2,
                        &eps, &domg12, numit < maxit1, &dv, Ca);
        /* 2 * tol0 is approximately 1 ulp for a number in [0, pi]. */
        /* Reversed test to allow escape with NaNs */
        if (tripb || !(fabs(v) >= (tripn ? 8 : 1) * tol0)) break;
        /* Update bracketing values */
        if (v > 0 && (numit > maxit1 || calp1/salp1 > calp1b/salp1b))
          { salp1b = salp1; calp1b = calp1; }
        else if (v < 0 && (numit > maxit1 || calp1/salp1 < calp1a/salp1a))
          { salp1a = salp1; calp1a = calp1; }
        if (numit < maxit1 && dv > 0) {
          real
            dalp1 = -v/dv;
          real
            sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
            nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
          if (nsalp1 > 0 && fabs(dalp1) < pi) {
            calp1 = calp1 * cdalp1 - salp1 * sdalp1;
            salp1 = nsalp1;
            norm2(&salp1, &calp1);
            /* In some regimes we don't get quadratic convergence because
             * slope -> 0.  So use convergence conditions based on epsilon
             * instead of sqrt(epsilon). */
            tripn = fabs(v) <= 16 * tol0;
            continue;
          }
        }
        /* Either dv was not positive or updated value was outside legal
         * range.  Use the midpoint of the bracket as the next estimate.
         * This mechanism is not needed for the WGS84 ellipsoid, but it does
         * catch problems with more eccentric ellipsoids.  Its efficacy is
         * such for the WGS84 test set with the starting guess set to alp1 =
         * 90deg:
         * the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
         * WGS84 and random input: mean = 4.74, sd = 0.99 */
        salp1 = (salp1a + salp1b)/2;
        calp1 = (calp1a + calp1b)/2;
        norm2(&salp1, &calp1);
        tripn = FALSE;
        tripb = (fabs(salp1a - salp1) + (calp1a - calp1) < tolb ||
                 fabs(salp1 - salp1b) + (calp1 - calp1b) < tolb);
      }
      Lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
              cbet1, cbet2, &s12x, &m12x, 0,
              outmask & GEOD_GEODESICSCALE ? &M12 : 0,
              outmask & GEOD_GEODESICSCALE ? &M21 : 0, Ca);
      m12x *= g->b;
      s12x *= g->b;
      a12 = sig12 / degree;
      if (outmask & GEOD_AREA) {
        /* omg12 = lam12 - domg12 */
        real sdomg12 = sin(domg12), cdomg12 = cos(domg12);
        somg12 = slam12 * cdomg12 - clam12 * sdomg12;
        comg12 = clam12 * cdomg12 + slam12 * sdomg12;
      }
    }
  }

  if (outmask & GEOD_DISTANCE)
    s12 = 0 + s12x;             /* Convert -0 to 0 */

  if (outmask & GEOD_REDUCEDLENGTH)
    m12 = 0 + m12x;             /* Convert -0 to 0 */

  if (outmask & GEOD_AREA) {
    real
      /* From Lambda12: sin(alp1) * cos(bet1) = sin(alp0) */
      salp0 = salp1 * cbet1,
      calp0 = hypotx(calp1, salp1 * sbet1); /* calp0 > 0 */
    real alp12;
    if (calp0 != 0 && salp0 != 0) {
      real
        /* From Lambda12: tan(bet) = tan(sig) * cos(alp) */
        ssig1 = sbet1, csig1 = calp1 * cbet1,
        ssig2 = sbet2, csig2 = calp2 * cbet2,
        k2 = sq(calp0) * g->ep2,
        eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2),
        /* Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0). */
        A4 = sq(g->a) * calp0 * salp0 * g->e2;
      real B41, B42;
      norm2(&ssig1, &csig1);
      norm2(&ssig2, &csig2);
      C4f(g, eps, Ca);
      B41 = SinCosSeries(FALSE, ssig1, csig1, Ca, nC4);
      B42 = SinCosSeries(FALSE, ssig2, csig2, Ca, nC4);
      S12 = A4 * (B42 - B41);
    } else
      /* Avoid problems with indeterminate sig1, sig2 on equator */
      S12 = 0;

    if (!meridian && somg12 > 1) {
      somg12 = sin(omg12); comg12 = cos(omg12);
    }

    if (!meridian &&
        /* omg12 < 3/4 * pi */
        comg12 > -(real)(0.7071) &&     /* Long difference not too big */
        sbet2 - sbet1 < (real)(1.75)) { /* Lat difference not too big */
      /* Use tan(Gamma/2) = tan(omg12/2)
       * * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
       * with tan(x/2) = sin(x)/(1+cos(x)) */
      real
        domg12 = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
      alp12 = 2 * atan2( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
                         domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) );
    } else {
      /* alp12 = alp2 - alp1, used in atan2 so no need to normalize */
      real
        salp12 = salp2 * calp1 - calp2 * salp1,
        calp12 = calp2 * calp1 + salp2 * salp1;
      /* The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
       * salp12 = -0 and alp12 = -180.  However this depends on the sign
       * being attached to 0 correctly.  The following ensures the correct
       * behavior. */
      if (salp12 == 0 && calp12 < 0) {
        salp12 = tiny * calp1;
        calp12 = -1;
      }
      alp12 = atan2(salp12, calp12);
    }
    S12 += g->c2 * alp12;
    S12 *= swapp * lonsign * latsign;
    /* Convert -0 to 0 */
    S12 += 0;
  }

  /* Convert calp, salp to azimuth accounting for lonsign, swapp, latsign. */
  if (swapp < 0) {
    swapx(&salp1, &salp2);
    swapx(&calp1, &calp2);
    if (outmask & GEOD_GEODESICSCALE)
      swapx(&M12, &M21);
  }

  salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
  salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

  if (psalp1) *psalp1 = salp1;
  if (pcalp1) *pcalp1 = calp1;
  if (psalp2) *psalp2 = salp2;
  if (pcalp2) *pcalp2 = calp2;

  if (outmask & GEOD_DISTANCE)
    *ps12 = s12;
  if (outmask & GEOD_REDUCEDLENGTH)
    *pm12 = m12;
  if (outmask & GEOD_GEODESICSCALE) {
    if (pM12) *pM12 = M12;
    if (pM21) *pM21 = M21;
  }
  if (outmask & GEOD_AREA)
    *pS12 = S12;

  /* Returned value in [0, 180] */
  return a12;
}

real geod_geninverse(const struct geod_geodesic* g,
                     real lat1, real lon1, real lat2, real lon2,
                     real* ps12, real* pazi1, real* pazi2,
                     real* pm12, real* pM12, real* pM21, real* pS12) {
  real salp1, calp1, salp2, calp2,
    a12 = geod_geninverse_int(g, lat1, lon1, lat2, lon2, ps12,
                              &salp1, &calp1, &salp2, &calp2,
                              pm12, pM12, pM21, pS12);
  if (pazi1) *pazi1 = atan2dx(salp1, calp1);
  if (pazi2) *pazi2 = atan2dx(salp2, calp2);
  return a12;
}

void geod_inverse(const struct geod_geodesic* g,
                  real lat1, real lon1, real lat2, real lon2,
                  real* ps12, real* pazi1, real* pazi2) {
  geod_geninverse(g, lat1, lon1, lat2, lon2, ps12, pazi1, pazi2, 0, 0, 0, 0);
}

real SinCosSeries(boolx sinp, real sinx, real cosx, const real c[], int n) {
  /* Evaluate
   * y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
   *            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
   * using Clenshaw summation.  N.B. c[0] is unused for sin series
   * Approx operation count = (n + 5) mult and (2 * n + 2) add */
  real ar, y0, y1;
  c += (n + sinp);              /* Point to one beyond last element */
  ar = 2 * (cosx - sinx) * (cosx + sinx); /* 2 * cos(2 * x) */
  y0 = n & 1 ? *--c : 0; y1 = 0;          /* accumulators for sum */
  /* Now n is even */
  n /= 2;
  while (n--) {
    /* Unroll loop x 2, so accumulators return to their original role */
    y1 = ar * y0 - y1 + *--c;
    y0 = ar * y1 - y0 + *--c;
  }
  return sinp
    ? 2 * sinx * cosx * y0      /* sin(2 * x) * y0 */
    : cosx * (y0 - y1);         /* cos(x) * (y0 - y1) */
}

void Lengths(const struct geod_geodesic* g,
             real eps, real sig12,
             real ssig1, real csig1, real dn1,
             real ssig2, real csig2, real dn2,
             real cbet1, real cbet2,
             real* ps12b, real* pm12b, real* pm0,
             real* pM12, real* pM21,
             /* Scratch area of the right size */
             real Ca[]) {
  real m0 = 0, J12 = 0, A1 = 0, A2 = 0;
  real Cb[nC];

  /* Return m12b = (reduced length)/b; also calculate s12b = distance/b,
   * and m0 = coefficient of secular term in expression for reduced length. */
  boolx redlp = pm12b || pm0 || pM12 || pM21;
  if (ps12b || redlp) {
    A1 = A1m1f(eps);
    C1f(eps, Ca);
    if (redlp) {
      A2 = A2m1f(eps);
      C2f(eps, Cb);
      m0 = A1 - A2;
      A2 = 1 + A2;
    }
    A1 = 1 + A1;
  }
  if (ps12b) {
    real B1 = SinCosSeries(TRUE, ssig2, csig2, Ca, nC1) -
      SinCosSeries(TRUE, ssig1, csig1, Ca, nC1);
    /* Missing a factor of b */
    *ps12b = A1 * (sig12 + B1);
    if (redlp) {
      real B2 = SinCosSeries(TRUE, ssig2, csig2, Cb, nC2) -
        SinCosSeries(TRUE, ssig1, csig1, Cb, nC2);
      J12 = m0 * sig12 + (A1 * B1 - A2 * B2);
    }
  } else if (redlp) {
    /* Assume here that nC1 >= nC2 */
    int l;
    for (l = 1; l <= nC2; ++l)
      Cb[l] = A1 * Ca[l] - A2 * Cb[l];
    J12 = m0 * sig12 + (SinCosSeries(TRUE, ssig2, csig2, Cb, nC2) -
                        SinCosSeries(TRUE, ssig1, csig1, Cb, nC2));
  }
  if (pm0) *pm0 = m0;
  if (pm12b)
    /* Missing a factor of b.
     * Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
     * accurate cancellation in the case of coincident points. */
    *pm12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) -
      csig1 * csig2 * J12;
  if (pM12 || pM21) {
    real csig12 = csig1 * csig2 + ssig1 * ssig2;
    real t = g->ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
    if (pM12)
      *pM12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
    if (pM21)
      *pM21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
  }
}

real Astroid(real x, real y) {
  /* Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
   * This solution is adapted from Geocentric::Reverse. */
  real k;
  real
    p = sq(x),
    q = sq(y),
    r = (p + q - 1) / 6;
  if ( !(q == 0 && r <= 0) ) {
    real
      /* Avoid possible division by zero when r = 0 by multiplying equations
       * for s and t by r^3 and r, resp. */
      S = p * q / 4,            /* S = r^3 * s */
      r2 = sq(r),
      r3 = r * r2,
      /* The discriminant of the quadratic equation for T3.  This is zero on
       * the evolute curve p^(1/3)+q^(1/3) = 1 */
      disc = S * (S + 2 * r3);
    real u = r;
    real v, uv, w;
    if (disc >= 0) {
      real T3 = S + r3, T;
      /* Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
       * of precision due to cancellation.  The result is unchanged because
       * of the way the T is used in definition of u. */
      T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); /* T3 = (r * t)^3 */
      /* N.B. cbrtx always returns the real root.  cbrtx(-8) = -2. */
      T = cbrtx(T3);            /* T = r * t */
      /* T can be zero; but then r2 / T -> 0. */
      u += T + (T != 0 ? r2 / T : 0);
    } else {
      /* T is complex, but the way u is defined the result is real. */
      real ang = atan2(sqrt(-disc), -(S + r3));
      /* There are three possible cube roots.  We choose the root which
       * avoids cancellation.  Note that disc < 0 implies that r < 0. */
      u += 2 * r * cos(ang / 3);
    }
    v = sqrt(sq(u) + q);              /* guaranteed positive */
    /* Avoid loss of accuracy when u < 0. */
    uv = u < 0 ? q / (v - u) : u + v; /* u+v, guaranteed positive */
    w = (uv - q) / (2 * v);           /* positive? */
    /* Rearrange expression for k to avoid loss of accuracy due to
     * subtraction.  Division by 0 not possible because uv > 0, w >= 0. */
    k = uv / (sqrt(uv + sq(w)) + w);   /* guaranteed positive */
  } else {               /* q == 0 && r <= 0 */
    /* y = 0 with |x| <= 1.  Handle this case directly.
     * for y small, positive root is k = abs(y)/sqrt(1-x^2) */
    k = 0;
  }
  return k;
}

real InverseStart(const struct geod_geodesic* g,
                  real sbet1, real cbet1, real dn1,
                  real sbet2, real cbet2, real dn2,
                  real lam12, real slam12, real clam12,
                  real* psalp1, real* pcalp1,
                  /* Only updated if return val >= 0 */
                  real* psalp2, real* pcalp2,
                  /* Only updated for short lines */
                  real* pdnm,
                  /* Scratch area of the right size */
                  real Ca[]) {
  real salp1 = 0, calp1 = 0, salp2 = 0, calp2 = 0, dnm = 0;

  /* Return a starting point for Newton's method in salp1 and calp1 (function
   * value is -1).  If Newton's method doesn't need to be used, return also
   * salp2 and calp2 and function value is sig12. */
  real
    sig12 = -1,               /* Return value */
    /* bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0] */
    sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
    cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
  real sbet12a;
  boolx shortline = cbet12 >= 0 && sbet12 < (real)(0.5) &&
    cbet2 * lam12 < (real)(0.5);
  real somg12, comg12, ssig12, csig12;
#if defined(__GNUC__) && __GNUC__ == 4 &&       \
  (__GNUC_MINOR__ < 6 || defined(__MINGW32__))
  /* Volatile declaration needed to fix inverse cases
   * 88.202499451857 0 -88.202499451857 179.981022032992859592
   * 89.262080389218 0 -89.262080389218 179.992207982775375662
   * 89.333123580033 0 -89.333123580032997687 179.99295812360148422
   * which otherwise fail with g++ 4.4.4 x86 -O3 (Linux)
   * and g++ 4.4.0 (mingw) and g++ 4.6.1 (tdm mingw). */
  {
    volatile real xx1 = sbet2 * cbet1;
    volatile real xx2 = cbet2 * sbet1;
    sbet12a = xx1 + xx2;
  }
#else
  sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
#endif
  if (shortline) {
    real sbetm2 = sq(sbet1 + sbet2), omg12;
    /* sin((bet1+bet2)/2)^2
     * =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2) */
    sbetm2 /= sbetm2 + sq(cbet1 + cbet2);
    dnm = sqrt(1 + g->ep2 * sbetm2);
    omg12 = lam12 / (g->f1 * dnm);
    somg12 = sin(omg12); comg12 = cos(omg12);
  } else {
    somg12 = slam12; comg12 = clam12;
  }

  salp1 = cbet2 * somg12;
  calp1 = comg12 >= 0 ?
    sbet12 + cbet2 * sbet1 * sq(somg12) / (1 + comg12) :
    sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);

  ssig12 = hypotx(salp1, calp1);
  csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

  if (shortline && ssig12 < g->etol2) {
    /* really short lines */
    salp2 = cbet1 * somg12;
    calp2 = sbet12 - cbet1 * sbet2 *
      (comg12 >= 0 ? sq(somg12) / (1 + comg12) : 1 - comg12);
    norm2(&salp2, &calp2);
    /* Set return value */
    sig12 = atan2(ssig12, csig12);
  } else if (fabs(g->n) > (real)(0.1) || /* No astroid calc if too eccentric */
             csig12 >= 0 ||
             ssig12 >= 6 * fabs(g->n) * pi * sq(cbet1)) {
    /* Nothing to do, zeroth order spherical approximation is OK */
  } else {
    /* Scale lam12 and bet2 to x, y coordinate system where antipodal point
     * is at origin and singular point is at y = 0, x = -1. */
    real y, lamscale, betscale;
    /* Volatile declaration needed to fix inverse case
     * 56.320923501171 0 -56.320923501171 179.664747671772880215
     * which otherwise fails with g++ 4.4.4 x86 -O3 */
    volatile real x;
    real lam12x = atan2(-slam12, -clam12); /* lam12 - pi */
    if (g->f >= 0) {            /* In fact f == 0 does not get here */
      /* x = dlong, y = dlat */
      {
        real
          k2 = sq(sbet1) * g->ep2,
          eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
        lamscale = g->f * cbet1 * A3f(g, eps) * pi;
      }
      betscale = lamscale * cbet1;

      x = lam12x / lamscale;
      y = sbet12a / betscale;
    } else {                    /* f < 0 */
      /* x = dlat, y = dlong */
      real
        cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
        bet12a = atan2(sbet12a, cbet12a);
      real m12b, m0;
      /* In the case of lon12 = 180, this repeats a calculation made in
       * Inverse. */
      Lengths(g, g->n, pi + bet12a,
              sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
              cbet1, cbet2, 0, &m12b, &m0, 0, 0, Ca);
      x = -1 + m12b / (cbet1 * cbet2 * m0 * pi);
      betscale = x < -(real)(0.01) ? sbet12a / x :
        -g->f * sq(cbet1) * pi;
      lamscale = betscale / cbet1;
      y = lam12x / lamscale;
    }

    if (y > -tol1 && x > -1 - xthresh) {
      /* strip near cut */
      if (g->f >= 0) {
        salp1 = minx((real)(1), -(real)(x)); calp1 = - sqrt(1 - sq(salp1));
      } else {
        calp1 = maxx((real)(x > -tol1 ? 0 : -1), (real)(x));
        salp1 = sqrt(1 - sq(calp1));
      }
    } else {
      /* Estimate alp1, by solving the astroid problem.
       *
       * Could estimate alpha1 = theta + pi/2, directly, i.e.,
       *   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
       *   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
       *
       * However, it's better to estimate omg12 from astroid and use
       * spherical formula to compute alp1.  This reduces the mean number of
       * Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
       * (min 0 max 5).  The changes in the number of iterations are as
       * follows:
       *
       * change percent
       *    1       5
       *    0      78
       *   -1      16
       *   -2       0.6
       *   -3       0.04
       *   -4       0.002
       *
       * The histogram of iterations is (m = number of iterations estimating
       * alp1 directly, n = number of iterations estimating via omg12, total
       * number of trials = 148605):
       *
       *  iter    m      n
       *    0   148    186
       *    1 13046  13845
       *    2 93315 102225
       *    3 36189  32341
       *    4  5396      7
       *    5   455      1
       *    6    56      0
       *
       * Because omg12 is near pi, estimate work with omg12a = pi - omg12 */
      real k = Astroid(x, y);
      real
        omg12a = lamscale * ( g->f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k );
      somg12 = sin(omg12a); comg12 = -cos(omg12a);
      /* Update spherical estimate of alp1 using omg12 instead of lam12 */
      salp1 = cbet2 * somg12;
      calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
    }
  }
  /* Sanity check on starting guess.  Backwards check allows NaN through. */
  if (!(salp1 <= 0))
    norm2(&salp1, &calp1);
  else {
    salp1 = 1; calp1 = 0;
  }

  *psalp1 = salp1;
  *pcalp1 = calp1;
  if (shortline)
    *pdnm = dnm;
  if (sig12 >= 0) {
    *psalp2 = salp2;
    *pcalp2 = calp2;
  }
  return sig12;
}

real Lambda12(const struct geod_geodesic* g,
              real sbet1, real cbet1, real dn1,
              real sbet2, real cbet2, real dn2,
              real salp1, real calp1,
              real slam120, real clam120,
              real* psalp2, real* pcalp2,
              real* psig12,
              real* pssig1, real* pcsig1,
              real* pssig2, real* pcsig2,
              real* peps,
              real* pdomg12,
              boolx diffp, real* pdlam12,
              /* Scratch area of the right size */
              real Ca[]) {
  real salp2 = 0, calp2 = 0, sig12 = 0,
    ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0,
    domg12 = 0, dlam12 = 0;
  real salp0, calp0;
  real somg1, comg1, somg2, comg2, somg12, comg12, lam12;
  real B312, eta, k2;

  if (sbet1 == 0 && calp1 == 0)
    /* Break degeneracy of equatorial line.  This case has already been
     * handled. */
    calp1 = -tiny;

  /* sin(alp1) * cos(bet1) = sin(alp0) */
  salp0 = salp1 * cbet1;
  calp0 = hypotx(calp1, salp1 * sbet1); /* calp0 > 0 */

  /* tan(bet1) = tan(sig1) * cos(alp1)
   * tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1) */
  ssig1 = sbet1; somg1 = salp0 * sbet1;
  csig1 = comg1 = calp1 * cbet1;
  norm2(&ssig1, &csig1);
  /* norm2(&somg1, &comg1); -- don't need to normalize! */

  /* Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
   * about this case, since this can yield singularities in the Newton
   * iteration.
   * sin(alp2) * cos(bet2) = sin(alp0) */
  salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
  /* calp2 = sqrt(1 - sq(salp2))
   *       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
   * and subst for calp0 and rearrange to give (choose positive sqrt
   * to give alp2 in [0, pi/2]). */
  calp2 = cbet2 != cbet1 || fabs(sbet2) != -sbet1 ?
    sqrt(sq(calp1 * cbet1) +
         (cbet1 < -sbet1 ?
          (cbet2 - cbet1) * (cbet1 + cbet2) :
          (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
    fabs(calp1);
  /* tan(bet2) = tan(sig2) * cos(alp2)
   * tan(omg2) = sin(alp0) * tan(sig2). */
  ssig2 = sbet2; somg2 = salp0 * sbet2;
  csig2 = comg2 = calp2 * cbet2;
  norm2(&ssig2, &csig2);
  /* norm2(&somg2, &comg2); -- don't need to normalize! */

  /* sig12 = sig2 - sig1, limit to [0, pi] */
  sig12 = atan2(maxx((real)(0), csig1 * ssig2 - ssig1 * csig2),
                                csig1 * csig2 + ssig1 * ssig2);

  /* omg12 = omg2 - omg1, limit to [0, pi] */
  somg12 = maxx((real)(0), comg1 * somg2 - somg1 * comg2);
  comg12 =                 comg1 * comg2 + somg1 * somg2;
  /* eta = omg12 - lam120 */
  eta = atan2(somg12 * clam120 - comg12 * slam120,
              comg12 * clam120 + somg12 * slam120);
  k2 = sq(calp0) * g->ep2;
  eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
  C3f(g, eps, Ca);
  B312 = (SinCosSeries(TRUE, ssig2, csig2, Ca, nC3-1) -
          SinCosSeries(TRUE, ssig1, csig1, Ca, nC3-1));
  domg12 = -g->f * A3f(g, eps) * salp0 * (sig12 + B312);
  lam12 = eta + domg12;

  if (diffp) {
    if (calp2 == 0)
      dlam12 = - 2 * g->f1 * dn1 / sbet1;
    else {
      Lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
              cbet1, cbet2, 0, &dlam12, 0, 0, 0, Ca);
      dlam12 *= g->f1 / (calp2 * cbet2);
    }
  }

  *psalp2 = salp2;
  *pcalp2 = calp2;
  *psig12 = sig12;
  *pssig1 = ssig1;
  *pcsig1 = csig1;
  *pssig2 = ssig2;
  *pcsig2 = csig2;
  *peps = eps;
  *pdomg12 = domg12;
  if (diffp)
    *pdlam12 = dlam12;

  return lam12;
}

real A3f(const struct geod_geodesic* g, real eps) {
  /* Evaluate A3 */
  return polyval(nA3 - 1, g->A3x, eps);
}

void C3f(const struct geod_geodesic* g, real eps, real c[]) {
  /* Evaluate C3 coeffs
   * Elements c[1] through c[nC3 - 1] are set */
  real mult = 1;
  int o = 0, l;
  for (l = 1; l < nC3; ++l) {   /* l is index of C3[l] */
    int m = nC3 - l - 1;        /* order of polynomial in eps */
    mult *= eps;
    c[l] = mult * polyval(m, g->C3x + o, eps);
    o += m + 1;
  }
}

void C4f(const struct geod_geodesic* g, real eps, real c[]) {
  /* Evaluate C4 coeffs
   * Elements c[0] through c[nC4 - 1] are set */
  real mult = 1;
  int o = 0, l;
  for (l = 0; l < nC4; ++l) {   /* l is index of C4[l] */
    int m = nC4 - l - 1;        /* order of polynomial in eps */
    c[l] = mult * polyval(m, g->C4x + o, eps);
    o += m + 1;
    mult *= eps;
  }
}

/* The scale factor A1-1 = mean value of (d/dsigma)I1 - 1 */
real A1m1f(real eps)  {
  static const real coeff[] = {
    /* (1-eps)*A1-1, polynomial in eps2 of order 3 */
    1, 4, 64, 0, 256,
  };
  int m = nA1/2;
  real t = polyval(m, coeff, sq(eps)) / coeff[m + 1];
  return (t + eps) / (1 - eps);
}

/* The coefficients C1[l] in the Fourier expansion of B1 */
void C1f(real eps, real c[])  {
  static const real coeff[] = {
    /* C1[1]/eps^1, polynomial in eps2 of order 2 */
    -1, 6, -16, 32,
    /* C1[2]/eps^2, polynomial in eps2 of order 2 */
    -9, 64, -128, 2048,
    /* C1[3]/eps^3, polynomial in eps2 of order 1 */
    9, -16, 768,
    /* C1[4]/eps^4, polynomial in eps2 of order 1 */
    3, -5, 512,
    /* C1[5]/eps^5, polynomial in eps2 of order 0 */
    -7, 1280,
    /* C1[6]/eps^6, polynomial in eps2 of order 0 */
    -7, 2048,
  };
  real
    eps2 = sq(eps),
    d = eps;
  int o = 0, l;
  for (l = 1; l <= nC1; ++l) {  /* l is index of C1p[l] */
    int m = (nC1 - l) / 2;      /* order of polynomial in eps^2 */
    c[l] = d * polyval(m, coeff + o, eps2) / coeff[o + m + 1];
    o += m + 2;
    d *= eps;
  }
}

/* The coefficients C1p[l] in the Fourier expansion of B1p */
void C1pf(real eps, real c[])  {
  static const real coeff[] = {
    /* C1p[1]/eps^1, polynomial in eps2 of order 2 */
    205, -432, 768, 1536,
    /* C1p[2]/eps^2, polynomial in eps2 of order 2 */
    4005, -4736, 3840, 12288,
    /* C1p[3]/eps^3, polynomial in eps2 of order 1 */
    -225, 116, 384,
    /* C1p[4]/eps^4, polynomial in eps2 of order 1 */
    -7173, 2695, 7680,
    /* C1p[5]/eps^5, polynomial in eps2 of order 0 */
    3467, 7680,
    /* C1p[6]/eps^6, polynomial in eps2 of order 0 */
    38081, 61440,
  };
  real
    eps2 = sq(eps),
    d = eps;
  int o = 0, l;
  for (l = 1; l <= nC1p; ++l) { /* l is index of C1p[l] */
    int m = (nC1p - l) / 2;     /* order of polynomial in eps^2 */
    c[l] = d * polyval(m, coeff + o, eps2) / coeff[o + m + 1];
    o += m + 2;
    d *= eps;
  }
}

/* The scale factor A2-1 = mean value of (d/dsigma)I2 - 1 */
real A2m1f(real eps)  {
  static const real coeff[] = {
    /* (eps+1)*A2-1, polynomial in eps2 of order 3 */
    -11, -28, -192, 0, 256,
  };
  int m = nA2/2;
  real t = polyval(m, coeff, sq(eps)) / coeff[m + 1];
  return (t - eps) / (1 + eps);
}

/* The coefficients C2[l] in the Fourier expansion of B2 */
void C2f(real eps, real c[])  {
  static const real coeff[] = {
    /* C2[1]/eps^1, polynomial in eps2 of order 2 */
    1, 2, 16, 32,
    /* C2[2]/eps^2, polynomial in eps2 of order 2 */
    35, 64, 384, 2048,
    /* C2[3]/eps^3, polynomial in eps2 of order 1 */
    15, 80, 768,
    /* C2[4]/eps^4, polynomial in eps2 of order 1 */
    7, 35, 512,
    /* C2[5]/eps^5, polynomial in eps2 of order 0 */
    63, 1280,
    /* C2[6]/eps^6, polynomial in eps2 of order 0 */
    77, 2048,
  };
  real
    eps2 = sq(eps),
    d = eps;
  int o = 0, l;
  for (l = 1; l <= nC2; ++l) { /* l is index of C2[l] */
    int m = (nC2 - l) / 2;     /* order of polynomial in eps^2 */
    c[l] = d * polyval(m, coeff + o, eps2) / coeff[o + m + 1];
    o += m + 2;
    d *= eps;
  }
}

/* The scale factor A3 = mean value of (d/dsigma)I3 */
void A3coeff(struct geod_geodesic* g) {
  static const real coeff[] = {
    /* A3, coeff of eps^5, polynomial in n of order 0 */
    -3, 128,
    /* A3, coeff of eps^4, polynomial in n of order 1 */
    -2, -3, 64,
    /* A3, coeff of eps^3, polynomial in n of order 2 */
    -1, -3, -1, 16,
    /* A3, coeff of eps^2, polynomial in n of order 2 */
    3, -1, -2, 8,
    /* A3, coeff of eps^1, polynomial in n of order 1 */
    1, -1, 2,
    /* A3, coeff of eps^0, polynomial in n of order 0 */
    1, 1,
  };
  int o = 0, k = 0, j;
  for (j = nA3 - 1; j >= 0; --j) {             /* coeff of eps^j */
    int m = nA3 - j - 1 < j ? nA3 - j - 1 : j; /* order of polynomial in n */
    g->A3x[k++] = polyval(m, coeff + o, g->n) / coeff[o + m + 1];
    o += m + 2;
  }
}

/* The coefficients C3[l] in the Fourier expansion of B3 */
void C3coeff(struct geod_geodesic* g) {
  static const real coeff[] = {
    /* C3[1], coeff of eps^5, polynomial in n of order 0 */
    3, 128,
    /* C3[1], coeff of eps^4, polynomial in n of order 1 */
    2, 5, 128,
    /* C3[1], coeff of eps^3, polynomial in n of order 2 */
    -1, 3, 3, 64,
    /* C3[1], coeff of eps^2, polynomial in n of order 2 */
    -1, 0, 1, 8,
    /* C3[1], coeff of eps^1, polynomial in n of order 1 */
    -1, 1, 4,
    /* C3[2], coeff of eps^5, polynomial in n of order 0 */
    5, 256,
    /* C3[2], coeff of eps^4, polynomial in n of order 1 */
    1, 3, 128,
    /* C3[2], coeff of eps^3, polynomial in n of order 2 */
    -3, -2, 3, 64,
    /* C3[2], coeff of eps^2, polynomial in n of order 2 */
    1, -3, 2, 32,
    /* C3[3], coeff of eps^5, polynomial in n of order 0 */
    7, 512,
    /* C3[3], coeff of eps^4, polynomial in n of order 1 */
    -10, 9, 384,
    /* C3[3], coeff of eps^3, polynomial in n of order 2 */
    5, -9, 5, 192,
    /* C3[4], coeff of eps^5, polynomial in n of order 0 */
    7, 512,
    /* C3[4], coeff of eps^4, polynomial in n of order 1 */
    -14, 7, 512,
    /* C3[5], coeff of eps^5, polynomial in n of order 0 */
    21, 2560,
  };
  int o = 0, k = 0, l, j;
  for (l = 1; l < nC3; ++l) {                    /* l is index of C3[l] */
    for (j = nC3 - 1; j >= l; --j) {             /* coeff of eps^j */
      int m = nC3 - j - 1 < j ? nC3 - j - 1 : j; /* order of polynomial in n */
      g->C3x[k++] = polyval(m, coeff + o, g->n) / coeff[o + m + 1];
      o += m + 2;
    }
  }
}

/* The coefficients C4[l] in the Fourier expansion of I4 */
void C4coeff(struct geod_geodesic* g) {
  static const real coeff[] = {
    /* C4[0], coeff of eps^5, polynomial in n of order 0 */
    97, 15015,
    /* C4[0], coeff of eps^4, polynomial in n of order 1 */
    1088, 156, 45045,
    /* C4[0], coeff of eps^3, polynomial in n of order 2 */
    -224, -4784, 1573, 45045,
    /* C4[0], coeff of eps^2, polynomial in n of order 3 */
    -10656, 14144, -4576, -858, 45045,
    /* C4[0], coeff of eps^1, polynomial in n of order 4 */
    64, 624, -4576, 6864, -3003, 15015,
    /* C4[0], coeff of eps^0, polynomial in n of order 5 */
    100, 208, 572, 3432, -12012, 30030, 45045,
    /* C4[1], coeff of eps^5, polynomial in n of order 0 */
    1, 9009,
    /* C4[1], coeff of eps^4, polynomial in n of order 1 */
    -2944, 468, 135135,
    /* C4[1], coeff of eps^3, polynomial in n of order 2 */
    5792, 1040, -1287, 135135,
    /* C4[1], coeff of eps^2, polynomial in n of order 3 */
    5952, -11648, 9152, -2574, 135135,
    /* C4[1], coeff of eps^1, polynomial in n of order 4 */
    -64, -624, 4576, -6864, 3003, 135135,
    /* C4[2], coeff of eps^5, polynomial in n of order 0 */
    8, 10725,
    /* C4[2], coeff of eps^4, polynomial in n of order 1 */
    1856, -936, 225225,
    /* C4[2], coeff of eps^3, polynomial in n of order 2 */
    -8448, 4992, -1144, 225225,
    /* C4[2], coeff of eps^2, polynomial in n of order 3 */
    -1440, 4160, -4576, 1716, 225225,
    /* C4[3], coeff of eps^5, polynomial in n of order 0 */
    -136, 63063,
    /* C4[3], coeff of eps^4, polynomial in n of order 1 */
    1024, -208, 105105,
    /* C4[3], coeff of eps^3, polynomial in n of order 2 */
    3584, -3328, 1144, 315315,
    /* C4[4], coeff of eps^5, polynomial in n of order 0 */
    -128, 135135,
    /* C4[4], coeff of eps^4, polynomial in n of order 1 */
    -2560, 832, 405405,
    /* C4[5], coeff of eps^5, polynomial in n of order 0 */
    128, 99099,
  };
  int o = 0, k = 0, l, j;
  for (l = 0; l < nC4; ++l) {        /* l is index of C4[l] */
    for (j = nC4 - 1; j >= l; --j) { /* coeff of eps^j */
      int m = nC4 - j - 1;           /* order of polynomial in n */
      g->C4x[k++] = polyval(m, coeff + o, g->n) / coeff[o + m + 1];
      o += m + 2;
    }
  }
}

int transit(real lon1, real lon2) {
  real lon12;
  /* Return 1 or -1 if crossing prime meridian in east or west direction.
   * Otherwise return zero. */
  /* Compute lon12 the same way as Geodesic::Inverse. */
  lon1 = AngNormalize(lon1);
  lon2 = AngNormalize(lon2);
  lon12 = AngDiff(lon1, lon2, 0);
  return lon1 <= 0 && lon2 > 0 && lon12 > 0 ? 1 :
    (lon2 <= 0 && lon1 > 0 && lon12 < 0 ? -1 : 0);
}

int transitdirect(real lon1, real lon2) {
#if HAVE_C99_MATH
  lon1 = remainder(lon1, (real)(720));
  lon2 = remainder(lon2, (real)(720));
  return ( (lon2 >= 0 && lon2 < 360 ? 0 : 1) -
           (lon1 >= 0 && lon1 < 360 ? 0 : 1) );
#else
  lon1 = fmod(lon1, (real)(720));
  lon2 = fmod(lon2, (real)(720));
  return ( ((lon2 >= 0 && lon2 < 360) || lon2 < -360 ? 0 : 1) -
           ((lon1 >= 0 && lon1 < 360) || lon1 < -360 ? 0 : 1) );
#endif
}

void accini(real s[]) {
  /* Initialize an accumulator; this is an array with two elements. */
  s[0] = s[1] = 0;
}

void acccopy(const real s[], real t[]) {
  /* Copy an accumulator; t = s. */
  t[0] = s[0]; t[1] = s[1];
}

void accadd(real s[], real y) {
  /* Add y to an accumulator. */
  real u, z = sumx(y, s[1], &u);
  s[0] = sumx(z, s[0], &s[1]);
  if (s[0] == 0)
    s[0] = u;
  else
    s[1] = s[1] + u;
}

real accsum(const real s[], real y) {
  /* Return accumulator + y (but don't add to accumulator). */
  real t[2];
  acccopy(s, t);
  accadd(t, y);
  return t[0];
}

void accneg(real s[]) {
  /* Negate an accumulator. */
  s[0] = -s[0]; s[1] = -s[1];
}

/** @endcond */
