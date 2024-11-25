#ifndef DISTS_SEQ_H
#define DISTS_SEQ_H

#include <R.h>
#include <Rinternals.h>

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_seq (SEXP x_);
SEXP R_vincenty_seq (SEXP x_);
SEXP R_cheap_seq (SEXP x_);
SEXP R_geodesic_seq (SEXP x_);

#endif /* DISTS_SEQ_H */
