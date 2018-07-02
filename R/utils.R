#' geodist_benchmark
#'
#' Benchmark errors for different geodist measures
#'
#' @param lat Central latgitude where errors should be measured
#' @param d Distance in metres over which errors should be measured
#' @param n Number of random values used to generate estimates
#' @export
#' @examples
#' geodist_benchmark (0, 1, 100)
geodist_benchmark <- function (lat = 0, d = 1, n = 1e2)
{
    lon <- 0
    dist_methods <- c ("geodesic", "haversine", "vincenty", "cheap")
    delta <- get_delta (lon, lat, d)
    x <- cbind ((lon - delta / 2) + delta * runif (n),
                (lat - delta / 2) + delta * runif (n))
    colnames (x) <- c ("x", "y")
    d <- lapply (dist_methods, function (i) {
                     res <- geodist (x, measure = i)
                     res [upper.tri (res)]
                })
    indx <- which (!grepl ("geodesic", dist_methods))
    dabs <- unlist (lapply (d [indx], function (i) mean (abs (i - d [[1]]))))
    drel <- unlist (lapply (d [indx], function (i)
                            mean (abs (i - d [[1]]) / d [[1]])))
    res <- rbind (dabs, drel)
    rownames (res) <- c ("absolute", "relative")
    colnames (res) <- dist_methods [indx]
    return (res)
}

#' get_delta
#'
#' Invert the inverse method to find the rough range of lon/lat corresponding to
#' a given distance. This is only used for geodist_benchmark, so the
#' computational inefficiency is no real issue here
#' @noRd
get_delta <- function (lon = 0, lat = 0, d = 1)
{
    f <- function (delta, lon, lat, d) {
        x <- c (lon - delta / 2, lat - delta / 2)
        y <- c (lon + delta / 2, lat + delta / 2)
        names (x) <- names (y) <- c ("x", "y")
        abs ((as.numeric (geodist (x, y, measure = "haversine")) - d))
    }
    optim (d / 1e5, f, lon = lon, lat = lat, d = d, method = "Brent",
           lower = 0, upper = d / 1e5)$par
}

#' find_xy_cols
#'
#' Find the lon and lat cols of a rectangular object
#' @param obj Rectangular object
#' @return Vector of two column indices of longitude and latitude
#' @noRd
find_xy_cols <- function (obj)
{
    nms <- names (obj)
    if (is.null (nms))
        nms <- colnames (obj)

    if (!is.null (nms))
    {
        ix <- which (grepl ("x", nms, ignore.case = TRUE) |
                     grepl ("lon", nms, ignore.case = TRUE))
        if (length (ix) > 0)
            ix <- which (grepl ("^x", nms, ignore.case = TRUE) |
                         grepl ("lon", nms, ignore.case = TRUE))
        iy <- which (grepl ("y", nms, ignore.case = TRUE) |
                     grepl ("lat", nms, ignore.case = TRUE))
        if (length (iy) > 0)
            iy <- which (grepl ("^y", nms, ignore.case = TRUE) |
                         grepl ("lat", nms, ignore.case = TRUE))
    } else
    {
        message ("object has no named columns; assuming order is lon then lat")
        ix <- 1
        iy <- 2
    }
    c (ix, iy)
}

#' convert_to_matrix
#'
#' Convert a rectangular object to a matrix of two columns, lon and lat
#'
#' @param obj Rectagular object
#' @return Numeric matrix of two columns: lon and lat
#' @noRd
convert_to_matrix <- function (obj)
{
    xy_cols <- find_xy_cols (obj)
    if (is.vector (obj))
        obj <- matrix (obj, nrow = 1)
    if (is.numeric (obj))
        cbind (as.numeric (obj [, xy_cols [1]]),
               as.numeric (obj [, xy_cols [2]]))
    else
        cbind (as.numeric (obj [[xy_cols [1] ]]),
               as.numeric (obj [[xy_cols [2] ]]))
}
