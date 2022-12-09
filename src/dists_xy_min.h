#ifndef DISTS_XY_MIN_H
#define DISTS_XY_MIN_H

#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_xy_min (SEXP x_, SEXP y_);
SEXP R_vincenty_xy_min (SEXP x_, SEXP y_);
SEXP R_cheap_xy_min (SEXP x_, SEXP y_);
SEXP R_geodesic_xy_min (SEXP x_, SEXP y_);

#endif /* DISTS_XY_MIN_H */
