#ifndef DISTS_X_VEC_H
#define DISTS_X_VEC_H

#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_vec (SEXP x_, SEXP y_);
SEXP R_vincenty_vec (SEXP x_, SEXP y_);
SEXP R_cheap_vec (SEXP x_, SEXP y_);
SEXP R_geodesic_vec (SEXP x_, SEXP y_);

#endif /* DISTS_X_VEC_H */
