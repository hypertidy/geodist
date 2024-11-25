#ifndef DISTS_PAIRED_H
#define DISTS_PAIRED_H

#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_paired (SEXP x_, SEXP y_);
SEXP R_vincenty_paired (SEXP x_, SEXP y_);
SEXP R_cheap_paired (SEXP x_, SEXP y_);
SEXP R_geodesic_paired (SEXP x_, SEXP y_);

#endif /* DISTS_PAIRED_H */
