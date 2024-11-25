#ifndef DISTS_SEQ_VEC_H
#define DISTS_SEQ_VEC_H

#include <R.h>
#include <Rinternals.h>

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_seq_vec (SEXP x_, SEXP y_);
SEXP R_vincenty_seq_vec (SEXP x_, SEXP y_);
SEXP R_cheap_seq_vec (SEXP x_, SEXP y_);
SEXP R_geodesic_seq_vec (SEXP x_, SEXP y_);

#endif /* DISTS_SEQ_VEC_H */
