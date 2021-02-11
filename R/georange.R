#' georange
#'
#' Calculate range of distances (min-max) between all points in one or two
#' rectangular objects containing lon-lat coordinates.
#'
#' @inheritParams geodist
#' @return A named vector of two numeric values: minimum and maximum, giving the
#' respective distances in metres.
#'
#' @note \code{measure = "cheap"} denotes the mapbox cheap ruler
#' \url{https://github.com/mapbox/cheap-ruler-cpp}; \code{measure = "geodesic"}
#' denotes the very accurate geodesic methods given in Karney (2013)
#' "Algorithms for geodesics" J Geod 87:43-55, and as provided by the
#' code{sf::st_dist()} function.
#'
#' @export
#'
#' @examples
#' n <- 50
#' x <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
#' y <- cbind (-10 + 20 * runif (2 * n), -10 + 20 * runif (2 * n))
#' colnames (x) <- colnames (y) <- c ("x", "y")
#' # All of the following returns vector of two values: minimum and maximum:
#' d0 <- georange (x)
#' d1 <- georange (x, y)
#' d2 <- georange (x, sequential = TRUE)
#' d0_2 <- georange (x, measure = "geodesic") # nanometre-accurate version of d0
georange <- function (x, y, sequential = FALSE, measure = "cheap") {

    measures <- c ("haversine", "vincenty", "cheap", "geodesic")
    measure <- match.arg (tolower (measure), measures)
    x <- convert_to_matrix (x)
    if (!missing (y)) {

        if (sequential) {

            message ("Sequential distances calculated along values of 'x' only")
            georange_seq (x, measure)
        } else {

            y <- convert_to_matrix (y)
            georange_xy (x, y, measure)
            # t() because the src code loops over x then y, so y is the internal
            # loop
        }
    } else {

        if (sequential)
            georange_seq (x, measure)
        else
            georange_x (x, measure)
    }
}

georange_seq <- function (x, measure) {

    if (measure == "haversine")
        res <- .Call ("R_haversine_seq_range", as.vector (x))
    else if (measure == "vincenty")
        res <- .Call ("R_vincenty_seq_range", as.vector (x))
    else if (measure == "geodesic")
        res <- .Call ("R_geodesic_seq_range", as.vector (x))
    else
        res <- .Call ("R_cheap_seq_range", as.vector (x))

    names (res) <- c ("minimum", "maximum")

    return (res)
}

georange_x <- function (x, measure) {

    if (measure == "haversine")
        res <- .Call ("R_haversine_range", as.vector (x))
    else if (measure == "vincenty")
        res <- .Call ("R_vincenty_range", as.vector (x))
    else if (measure == "geodesic")
        res <- .Call ("R_geodesic_range", as.vector (x))
    else
        res <- .Call ("R_cheap_range", as.vector (x))

    names (res) <- c ("minimum", "maximum")

    return (res)
}

georange_xy <- function (x, y, measure) {

    if (measure == "haversine")
        res <- .Call ("R_haversine_xy_range", as.vector (x), as.vector (y))
    else if (measure == "vincenty")
        res <- .Call ("R_vincenty_xy_range", as.vector (x), as.vector (y))
    else if (measure == "geodesic")
        res <- .Call ("R_geodesic_xy_range", as.vector (x), as.vector (y))
    else if (measure == "cheap")
        res <- .Call ("R_cheap_xy_range", as.vector (x), as.vector (y))

    names (res) <- c ("minimum", "maximum")

    return (res)
}
