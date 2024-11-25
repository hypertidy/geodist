#ifndef DISTS_PAIRED_VEC_H
#define DISTS_PAIRED_VEC_H

#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_paired_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_);
SEXP R_vincenty_paired_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_);
SEXP R_cheap_paired_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_);
SEXP R_geodesic_paired_vec (SEXP x1_, SEXP y1_, SEXP x2_, SEXP y2_);

#endif /* DISTS_PAIRED_VEC_H */
