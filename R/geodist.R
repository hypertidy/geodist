#' geodist
#'
#' Convert one or two rectangular objects containing lon-lat coordinates into
#' vector or matrix of geodesic distances.
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
#' @return If only \code{x} passed and \code{sequential = FALSE}, a square
#' symmetric matrix containing distances between all items in \code{x}; If only
#' \code{x} passed and \code{sequential = TRUE}, a vector of sequential
#' distances between rows of \code{x}; otherwise if \code{y} is passed, a matrix
#' of \code{nrow(x)} rows and \code{nrow(y)} columns.
#'
#' @note \code{measure = "cheap"} denotes the mapbox cheap ruler
#' \url{https://github.com/mapbox/cheap-ruler-cpp}; \code{measure = "geodetic"}
#' denotes the very accurate geodetic methods given in Kearney (2013)
#' "Algorithms for geodesics" J Geod 87:43-55, and as provided by the 
#' code{sf::st_dist()} function.
#'
#' @export
#' @useDynLib geodist R_haversine R_vincenty R_cheap R_geodesic
#' @useDynLib geodist R_haversine_xy R_vincenty_xy R_cheap_xy R_geodesic_xy
#' @useDynLib geodist R_haversine_paired R_vincenty_paired R_cheap_paired R_geodesic_paired
#' @useDynLib geodist R_haversine_seq R_vincenty_seq R_cheap_seq R_geodesic_seq
#'
#' @examples
#' n <- 50
#' x <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
#' y <- cbind (-10 + 20 * runif (2 * n), -10 + 20 * runif (2 * n))
#' colnames (x) <- colnames (y) <- c ("x", "y")
#' d0 <- geodist (x) # A 50-by-50 matrix
#' d1 <- geodist (x, y) # A 50-by-100 matrix
#' d2 <- geodist (x, sequential = TRUE) # Vector of length 49
#' d2 <- geodist (x, sequential = TRUE, pad = TRUE) # Vector of length 50
#' d0_2 <- geodist (x, measure = "geodesic") # nanometre-accurate version of d0
geodist <- function (x, y, paired = FALSE,
                     sequential = FALSE, pad = FALSE, measure = "cheap")
{
    measures <- c ("haversine", "vincenty", "cheap", "geodesic")
    measure <- match.arg (tolower (measure), measures)
    x <- convert_to_matrix (x)
    if (!missing (y))
    {
        if (paired)
        {
            if (nrow (x) != nrow (y))
                stop ("x and y must have the same number of ",
                      "rows for paired distances")
            y <- convert_to_matrix (y)
            geodist_paired (x, y, measure)
        } else if (sequential)
        {
            message ("Sequential distances calculated along values of 'x' only")
            geodist_seq (x, measure, pad)
        } else
        {
            y <- convert_to_matrix (y)
            geodist_xy (x, y, measure)
            # t() because the src code loops over x then y, so y is the internal
            # loop
        }
    } else
    {
        if (sequential)
            geodist_seq (x, measure, pad)
        else
            geodist_x (x, measure)
    }
}

geodist_paired <- function (x, y, measure)
{
    if (measure == "haversine")
        res <- .Call("R_haversine_paired", as.vector (x), as.vector (y))
    else if (measure == "vincenty")
        res <- .Call("R_vincenty_paired", as.vector (x), as.vector (y))
    else if (measure == "geodesic")
        res <- .Call("R_geodesic_paired", as.vector (x), as.vector (y))
    else
        res <- .Call("R_cheap_paired", as.vector (x), as.vector (y))

    return (res)
}

geodist_seq <- function (x, measure, pad)
{
    if (measure == "haversine")
        res <- matrix (.Call("R_haversine_seq", as.vector (x)), nrow = nrow (x))
    else if (measure == "vincenty")
        res <- matrix (.Call("R_vincenty_seq", as.vector (x)), nrow = nrow (x))
    else if (measure == "geodesic")
        res <- matrix (.Call("R_geodesic_seq", as.vector (x)), nrow = nrow (x))
    else
        res <- matrix (.Call("R_cheap_seq", as.vector (x)), nrow = nrow (x))
    indx <- 1:length (res)
    if (!pad)
        indx <- 2:length (res)

    return (res [indx]) # implicitly converts to vector
}

geodist_x <- function (x, measure)
{
    if (measure == "haversine")
        matrix (.Call("R_haversine", as.vector (x)), nrow = nrow (x))
    else if (measure == "vincenty")
        matrix (.Call("R_vincenty", as.vector (x)), nrow = nrow (x))
    else if (measure == "geodesic")
        matrix (.Call("R_geodesic", as.vector (x)), nrow = nrow (x))
    else
        matrix (.Call("R_cheap", as.vector (x)), nrow = nrow (x))
}

geodist_xy <- function (x, y, measure)
{
    if (measure == "haversine")
        res <- .Call ("R_haversine_xy", as.vector (x), as.vector (y))
    else if (measure == "vincenty")
        res <- .Call ("R_vincenty_xy", as.vector (x), as.vector (y))
    else if (measure == "geodesic")
        res <- .Call("R_geodesic_xy", as.vector (x), as.vector (y))
    else if (measure == "cheap")
        res <- .Call ("R_cheap_xy", as.vector (x), as.vector (y))
    t (matrix (res, nrow = nrow (y)))
}
