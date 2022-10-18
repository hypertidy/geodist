#ifndef DISTS_X_H
#define DISTS_X_H

#include <R.h>
#include <Rinternals.h>

#include <stdio.h> 

#include "common.h"
#include "WSG84-defs.h"

SEXP R_haversine (SEXP x_);
SEXP R_vincenty (SEXP x_);
SEXP R_cheap (SEXP x_);
SEXP R_geodesic (SEXP x_);

#endif /* DISTS_X_H */
