#' geodist
#'
#' Convert one or two rectangular objects containing lon-lat coordinates into
#' matrix of geodesic distances.
#'
#' @param x Rectangular object (matrix, \code{data.frame}, \pkg{tibble},
#' whatever) containing longitude and latitude coordinates.
#' @param y Optional second object which, if passed, results in distances
#' calculated between each object in \code{x} and each in \code{y}.
#' @return If only \code{x} passed, a square matrix containing distances between
#' all items in \code{x}; otherwise if \code{y} is passed, the resultant matrix
#' has \code{nrow(x)} rows and \code{nrow(y)} columns.
#'
#' @export
#' @useDynLib geodist R_haversine
#' @useDynLib geodist R_haversine_xy
geodist <- function (x, y)
{
    x <- convert_to_matrix (x)
    if (!missing (y))
    {
        y <- convert_to_matrix (y)
        # t() because the src code loops over x then y, so y is the internal
        # loop
        t (matrix (.Call ("R_haversine_xy", as.vector (x), as.vector (y)),
                   nrow = nrow (y)))
    } else
    {
        matrix (.Call("R_haversine", as.vector (x)), nrow = nrow (x))
    }
}
