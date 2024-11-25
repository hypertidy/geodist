#' geodist
#'
#' Convert one or two rectangular objects containing lon-lat coordinates into
#' vector or matrix of geodesic distances in metres.
#'
#' @param x Rectangular object (matrix, \code{data.frame}, \pkg{tibble},
#' whatever) containing longitude and latitude coordinates.
#' @param y Optional second object which, if passed, results in distances
#' calculated between each object in \code{x} and each in \code{y}.
#' @param paired If \code{TRUE}, calculate paired distances between each entry
#' in \code{x} and \code{y}, returning a single vector.
#' @param sequential If \code{TRUE}, calculate (vector of) distances
#' sequentially along \code{x} (when no \code{y} is passed), otherwise calculate
#' matrix of pairwise distances between all points.
#' @param pad If \code{sequential = TRUE} values are padded with initial
#' \code{NA} to return \code{n} values for input with \code{n} rows, otherwise
#' return \code{n - 1} values.
#' @param measure One of "haversine" "vincenty", "geodesic", or "cheap"
#' specifying desired method of geodesic distance calculation; see Notes.
#' @param quiet If \code{FALSE}, check whether max of calculated distances
#' is greater than accuracy threshold and warn.
#' @return If only \code{x} passed and \code{sequential = FALSE}, a square
#' symmetric matrix containing distances between all items in \code{x}; If only
#' \code{x} passed and \code{sequential = TRUE}, a vector of sequential
#' distances between rows of \code{x}; otherwise if \code{y} is passed, a matrix
#' of \code{nrow(x)} rows and \code{nrow(y)} columns. All return values are
#' distances in metres.
#'
#' @note \code{measure = "cheap"} denotes the mapbox cheap ruler
#' \url{https://github.com/mapbox/cheap-ruler-cpp}; \code{measure = "geodesic"}
#' denotes the very accurate geodesic methods given in Karney (2013)
#' "Algorithms for geodesics" J Geod 87:43-55, and as provided by the
#' `st_dist()` function from the \pkg{sf} package.
#'
#' @export
#'
#' @examples
#' n <- 50
#' # Default "cheap" distance measure is only accurate for short distances:
#' x <- cbind (runif (n, -0.1, 0.1), runif (n, -0.1, 0.1))
#' y <- cbind (runif (2 * n, -0.1, 0.1), runif (2 * n, -0.1, 0.1))
#' colnames (x) <- colnames (y) <- c ("x", "y")
#' d0 <- geodist (x) # A 50-by-50 matrix
#' d1 <- geodist (x, y) # A 50-by-100 matrix
#' d2 <- geodist (x, sequential = TRUE) # Vector of length 49
#' d2 <- geodist (x, sequential = TRUE, pad = TRUE) # Vector of length 50
#' d0_2 <- geodist (x, measure = "geodesic") # nanometre-accurate version of d0
#'
#' # Input data can also be 'data.frame' objects:
#' xy <- data.frame (x = runif (n, -0.1, 0.1), y = runif (n, -0.1, 0.1))
#' d <- geodist (xy)
geodist <- function (x, y, paired = FALSE,
                     sequential = FALSE, pad = FALSE,
                     measure = "cheap", quiet = FALSE) {

    measures <- c ("haversine", "vincenty", "cheap", "geodesic")
    measure <- match.arg (tolower (measure), measures)

    x <- convert_to_matrix (x)

    if (!missing (y)) {

        if (paired) {

            if (nrow (x) != nrow (y)) {
                stop (
                    "x and y must have the same number of ",
                    "rows for paired distances"
                )
            }
            y <- convert_to_matrix (y)
            res <- geodist_paired (x, y, measure)
        } else if (sequential) {

            message ("Sequential distances calculated along values of 'x' only")
            res <- geodist_seq (x, measure, pad)
        } else {

            y <- convert_to_matrix (y)
            res <- geodist_xy (x, y, measure)
            # t() because the src code loops over x then y, so y is the internal
            # loop
        }
    } else {

        if (sequential) {
            res <- geodist_seq (x, measure, pad)
        } else {
            res <- geodist_x (x, measure)
        }
    }

    if (measure == "cheap" & !quiet) {
        check_max_d (res, measure)
    }

    return (res)
}

geodist_paired <- function (x, y, measure) {

    if (measure == "haversine") {
        .Call ("R_haversine_paired", as.vector (x), as.vector (y))
    } else if (measure == "vincenty") {
        .Call ("R_vincenty_paired", as.vector (x), as.vector (y))
    } else if (measure == "geodesic") {
        .Call ("R_geodesic_paired", as.vector (x), as.vector (y))
    } else {
        .Call ("R_cheap_paired", as.vector (x), as.vector (y))
    }
}

geodist_seq <- function (x, measure, pad) {

    if (measure == "haversine") {
        res <- matrix (.Call ("R_haversine_seq", as.vector (x)),
            nrow = nrow (x)
        )
    } else if (measure == "vincenty") {
        res <- matrix (.Call ("R_vincenty_seq", as.vector (x)), nrow = nrow (x))
    } else if (measure == "geodesic") {
        res <- matrix (.Call ("R_geodesic_seq", as.vector (x)), nrow = nrow (x))
    } else {
        res <- matrix (.Call ("R_cheap_seq", as.vector (x)), nrow = nrow (x))
    }

    index <- seq_along (res)
    if (!pad) {
        index <- index [-1]
    }

    return (res [index]) # implicitly converts to vector
}

geodist_x <- function (x, measure) {

    if (measure == "haversine") {
        matrix (.Call ("R_haversine", as.vector (x)), nrow = nrow (x))
    } else if (measure == "vincenty") {
        matrix (.Call ("R_vincenty", as.vector (x)), nrow = nrow (x))
    } else if (measure == "geodesic") {
        matrix (.Call ("R_geodesic", as.vector (x)), nrow = nrow (x))
    } else {
        matrix (.Call ("R_cheap", as.vector (x)), nrow = nrow (x))
    }
}

geodist_xy <- function (x, y, measure) {

    if (measure == "haversine") {
        res <- .Call ("R_haversine_xy", as.vector (x), as.vector (y))
    } else if (measure == "vincenty") {
        res <- .Call ("R_vincenty_xy", as.vector (x), as.vector (y))
    } else if (measure == "geodesic") {
        res <- .Call ("R_geodesic_xy", as.vector (x), as.vector (y))
    } else if (measure == "cheap") {
        res <- .Call ("R_cheap_xy", as.vector (x), as.vector (y))
    }

    t (matrix (res, nrow = nrow (y)))
}
