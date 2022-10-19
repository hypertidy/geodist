#ifndef RANGE_SEQ_H
#define RANGE_SEQ_H

#include <R.h>
#include <Rinternals.h>

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine_seq_range (SEXP x_);
SEXP R_vincenty_seq_range (SEXP x_);
SEXP R_cheap_seq_range (SEXP x_);
SEXP R_geodesic_seq_range (SEXP x_);

#endif /* RANGE_SEQ_H */
