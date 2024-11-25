#ifndef DISTS_XY_H
#define DISTS_XY_H

#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_xy (SEXP x_, SEXP y_);
SEXP R_vincenty_xy (SEXP x_, SEXP y_);
SEXP R_cheap_xy (SEXP x_, SEXP y_);
SEXP R_geodesic_xy (SEXP x_, SEXP y_);

#endif /* DISTS_XY_H */
