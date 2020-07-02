#' geodist_vec
#'
#' An alternative interface to the main \link{geodist} function that directly
#' accepts inputs as individual vectors of coordinates, rather than the matrix
#' or `data.frame` inputs of the main function. This interface is provided for
#' cases where computational efficiency is important, and will generally provide
#' faster results than the main function.
#'
#' @param x1 Numeric vector of longitude coordinates
#' @param y1 Numeric vector of latitude coordinates
#' @param x2 Optional second numeric vector of longitude coordinates
#' @param y2 Optional second numeric vector of latitude coordinates
#' @param paired If \code{TRUE}, calculate paired distances between each entry
#' in \code{(x1, y1)} and \code{(x2, y2)}, returning a single vector.
#' @param sequential If \code{TRUE}, calculate (vector of) distances
#' sequentially along \code{(x1, y1)} (when no \code{(x2, y2)} are passed),
#' otherwise calculate matrix of pairwise distances between all points.
#' @param pad If \code{sequential = TRUE} values are padded with initial
#' \code{NA} to return \code{n} values for inputs of lenght \code{n}, otherwise
#' return \code{n - 1} values.
#' @param measure One of "haversine" "vincenty", "geodesic", or "cheap"
#' specifying desired method of geodesic distance calculation; see Notes.
#' @return If only \code{(x1, y1)} are passed and \code{sequential = FALSE}, a
#' square symmetric matrix containing distances between all items in \code{(x1,
#' y1)}; If only \code{(x1, y1)} are passed and \code{sequential = TRUE}, a
#' vector of sequential distances between matching elements of \code{(x1, y1)};
#' otherwise if \code{(x2, y2)} are passed, a matrix of \code{lenght(x1) ==
#' length(y1)} rows and \code{length(x2) == length(y2)} columns.
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
#' # Default "cheap" distance measure is only accurate for short distances:
#' x1 <- -1 + 2 * runif (n, -0.1, 0.1)
#' y1 <- -1 + 2 * runif (n, -0.1, 0.1)
#' d0 <- geodist_vec (x1, y1) # A 50-by-50 matrix
#' d2 <- geodist_vec (x1, y1, sequential = TRUE) # Vector of length 49
#' d2 <- geodist_vec (x1, y1, sequential = TRUE, pad = TRUE) # Vector of length 50
#' x2 <- -10 + 20 * runif (2 * n, -0.1, 0.1)
#' y2 <- -10 + 20 * runif (2 * n, -0.1, 0.1)
#' d1 <- geodist_vec (x1, y1, x2, y2) # A 50-by-100 matrix
geodist_vec <- function (x1, y1, x2, y2, paired = FALSE,
                         sequential = FALSE, pad = FALSE, measure = "cheap")
{
    measures <- c ("haversine", "vincenty", "cheap", "geodesic")
    measure <- match.arg (tolower (measure), measures)

    check_vec_inputs (x1, y1, 1)
    
    if (!missing (x2))
    {
        check_vec_inputs (x2, y2, 2)
        if (paired)
        {
            res <- geodist_paired_vec (x1, y1, x2, y2, measure)
        } else if (sequential)
        {
            message ("Sequential distances calculated along values of 'x' only")
            res <- geodist_seq_vec (x1, y2, measure, pad)
        } else
        {
            res <- geodist_xy_vec (x1, y1, x2, y2, measure)
        }
    } else
    {
        if (sequential)
            res <- geodist_seq_vec (x1, y1, measure, pad)
        else
            res <- geodist_x_vec (x1, y1, measure)
    }

    if (measure == "cheap")
        check_max_d (res, measure)

    return (res)
}

check_vec_inputs <- function (x, y, n = 1)
{
    if (missing (x) | missing (y))
        stop (paste0 ("x", n, " and y", n, " must be provided"))
    if (!(is.vector (x) && is.vector (y)))
        stop ("geodist_vec only accepts vector inputs")
    if (length (x) != length (y))
        stop (paste0 ("x", n, " and y", n, " must have the same length"))
}

geodist_paired_vec <- function (x1, y1, x2, y2, measure)
{
    if (measure == "haversine")
        .Call ("R_haversine_paired_vec", x1, y1, x2, y2)
    else if (measure == "vincenty")
        .Call ("R_vincenty_paired_vec", x1, y1, x2, y2)
    else if (measure == "geodesic")
        .Call ("R_geodesic_paired_vec", x1, y1, x2, y2)
    else
        .Call ("R_cheap_paired_vec", x1, y1, x2, y2)
}

geodist_seq_vec <- function (x, y, measure, pad)
{
    if (measure == "haversine")
        res <- matrix (.Call ("R_haversine_seq_vec", x, y))
    else if (measure == "vincenty")
        res <- matrix (.Call ("R_vincenty_seq_vec", x, y))
    else if (measure == "geodesic")
        res <- matrix (.Call ("R_geodesic_seq_vec", x, y))
    else
        res <- matrix (.Call ("R_cheap_seq_vec", x, y))
    
    indx <- 1:length (res)
    if (!pad)
        indx <- 2:length (res)

    return (res [indx]) # implicitly converts to vector
}

geodist_x_vec <- function (x, y, measure)
{
    if (measure == "haversine")
        matrix (.Call ("R_haversine_vec", x, y), nrow = length (x))
    else if (measure == "vincenty")
        matrix (.Call ("R_vincenty_vec", x, y), nrow = length (x))
    else if (measure == "geodesic")
        matrix (.Call ("R_geodesic_vec", x, y), nrow = length (x))
    else
        matrix (.Call ("R_cheap_vec", x, y), nrow = length (x))
}

geodist_xy_vec <- function (x1, y1, x2, y2, measure)
{
    if (measure == "haversine")
        res <- .Call ("R_haversine_xy_vec", x1, y1, x2, y2)
    else if (measure == "vincenty")
        res <- .Call ("R_vincenty_xy_vec", x1, y1, x2, y2)
    else if (measure == "geodesic")
        res <- .Call ("R_geodesic_xy_vec", x1, y1, x2, y2)
    else if (measure == "cheap")
        res <- .Call ("R_cheap_xy_vec", x1, y1, x2, y2)
    t (matrix (res, nrow = length (x2)))
}
