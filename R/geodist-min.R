#' Minimal pairwise distances between two input matrices
#'
#' Convert two rectangular objects containing lon-lat coordinates into an index
#' vector of the elements in the second matrix corresponding to the minimal
#' distance to each element of the first matrix.
#'
#' @param x Rectangular object (matrix, \code{data.frame}, \pkg{tibble},
#' whatever) containing longitude and latitude coordinates.
#' @param y Second rectangular object to be search for minimal distance to each
#' row in the first object.
#' @param measure One of "haversine" "vincenty", "geodesic", or "cheap"
#' specifying desired method of geodesic distance calculation; see Notes.
#' @param quiet If \code{FALSE}, check whether max of calculated distances
#' is greater than accuracy threshold and warn.
#' @return A integer index vector indexing elements of 'y' corresponding to
#' minimal distances to each element of 'x'. The length of this vector is equal
#' to the number of rows in 'x'.
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
#' index <- geodist_min (x, y, measure = "Haversine")
#' # 'index' is a vector of 50 values indexing elements of `y` corresponding to
#' # minimal distances to each element of `x`. It is exactly the same as the
#' # value returned by these lines:
#' d0 <- geodist (x, y, measure = "Haversine")
#' index0 <- apply (d0, 1, which.min)
#' identical (index, index0)
geodist_min <- function (x, y, measure = "cheap", quiet = FALSE) {

    measures <- c ("haversine", "vincenty", "cheap", "geodesic")
    measure <- match.arg (tolower (measure), measures)

    x <- convert_to_matrix (x)
    y <- convert_to_matrix (y)

    if (measure == "haversine") {
        res <- .Call ("R_haversine_xy_min", as.vector (x), as.vector (y))
    } else if (measure == "vincenty") {
        res <- .Call ("R_vincenty_xy_min", as.vector (x), as.vector (y))
    } else if (measure == "geodesic") {
        res <- .Call ("R_geodesic_xy_min", as.vector (x), as.vector (y))
    } else {
        res <- .Call ("R_cheap_xy_min", as.vector (x), as.vector (y))
    }

    return (res)
}
