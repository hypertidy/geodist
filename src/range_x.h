#ifndef RANGE_X_H
#define RANGE_X_H

#include <R.h>
#include <Rinternals.h>

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_range (SEXP x_);
SEXP R_vincenty_range (SEXP x_);
SEXP R_cheap_range (SEXP x_);
SEXP R_geodesic_range (SEXP x_);

#endif /* RANGE_X_H */
