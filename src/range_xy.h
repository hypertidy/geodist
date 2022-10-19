#ifndef RANGE_XY_H
#define RANGE_XY_H

#include <R.h>
#include <Rinternals.h>

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_xy_range (SEXP x_, SEXP y_);
SEXP R_vincenty_xy_range (SEXP x_, SEXP y_);
SEXP R_cheap_xy_range (SEXP x_, SEXP y_);
SEXP R_geodesic_xy_range (SEXP x_, SEXP y_);

#endif /* RANGE_XY_H */
