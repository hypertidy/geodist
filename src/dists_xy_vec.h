#ifndef DISTS_XY_H
#define DISTS_XY_H

#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_xy_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_);
SEXP R_vincenty_xy_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_);
SEXP R_cheap_xy_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_);
SEXP R_geodesic_xy_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_);

#endif /* DISTS_XY_H */
