#' geodist_benchmark
#'
#' Benchmark errors for different geodist measures
#'
#' @param lat Central latitude where errors should be measured
#' @param d Distance in metres over which errors should be measured
#' @param n Number of random values used to generate estimates
#' @return A 'data.frame' with three columns respectively comparing the accuracy
#' of the [Haversine, Vincenty, cheap] metrics against geodesic measures in both
#' absolute and relative terms (as two rows of the table).
#' @export
#' @examples
#' geodist_benchmark (0.0, 1.0, 100L)
geodist_benchmark <- function (lat = 0.0, d = 1.0, n = 100L) {

    chk_is_num_len_1 (lat)
    chk_is_num_len_1 (d)
    chk_is_num_len_1 (n)

    if (n < 2) {
        stop ("Comparisons require at least n = 2 objects")
    }
    if (n > 1e3) {
        stop (
            "benchmarking compares n ^ 2 estimates, and there's ",
            "nothing to be gained by extending beyond a million comparisons"
        )
    }

    lon <- 0
    dist_methods <- c ("geodesic", "haversine", "vincenty", "cheap")
    delta <- get_delta (lon, lat, d)

    x <- cbind (
        (lon - delta / 2) + delta * runif (n),
        (lat - delta / 2) + delta * runif (n)
    )
    colnames (x) <- c ("x", "y")

    d <- lapply (dist_methods, function (i) {
        res <- geodist (x, measure = i)
        res [upper.tri (res)]
    })

    indx <- which (!grepl ("geodesic", dist_methods))
    dabs <- unlist (lapply (d [indx], function (i) mean (abs (i - d [[1]]))))
    drel <- unlist (lapply (d [indx], function (i) {
        mean (abs (i - d [[1]]) / d [[1]])
    }))

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
get_delta <- function (lon = 0, lat = 0, d = 1) {

    f <- function (delta, lon, lat, d) {
        x <- c (lon - delta / 2, lat - delta / 2)
        y <- c (lon + delta / 2, lat + delta / 2)
        names (x) <- names (y) <- c ("x", "y")
        abs ((as.numeric (geodist (x, y, measure = "haversine")) - d))
    }

    optim (d / 1e5, f,
        lon = lon, lat = lat, d = d, method = "Brent",
        lower = 0, upper = d / 1e5
    )$par
}

#' find_xy_cols
#'
#' Find the lon and lat cols of a rectangular object
#' @param obj Rectangular object
#' @return Vector of two column indices of longitude and latitude
#' @noRd
find_xy_cols <- function (obj) {

    nms <- names (obj)
    if (is.null (nms)) {
        nms <- colnames (obj)
    }

    if (!is.null (nms)) {

        ix <- match_xy_pattern (nms, value = "x")
        iy <- match_xy_pattern (nms, value = "y")
        if (length (ix) != 1 | length (iy) != 1) {
            stop (
                "Unable to determine longitude and latitude columns; ",
                "perhaps try re-naming columns."
            )
        }
    } else {

        message ("object has no named columns; assuming order is lon then lat")
        ix <- 1
        iy <- 2
    }
    c (ix, iy)
}

match_xy_pattern <- function (nms, value = "x") {

    value <- match.arg (value, c ("x", "y"))
    lonlat <- ifelse (value == "x", "lon", "lat")

    ptn <- paste0 ("^", value, "|", value, "$|^", lonlat, "|", lonlat, "$")
    i <- grep (ptn, nms, ignore.case = TRUE)

    if (length (i) > 1) {
        # exclude any with :alpha: before or after x/lon:
        ptn <- paste0 (
            "^", value, "[[:alpha:]]|",
            "[[:alpha:]]", value, "$|",
            "^", lonlat, "[[:alpha:]]|",
            "[[:alpha:]]", lonlat, "$"
        )
        i <- i [which (!seq_along (i) %in% grep (ptn, nms [i]))]
    }

    if (length (i) == 0) {
        # search for match to initial part only
        ptn <- paste0 ("^", value, "|^", lonlat)
        i <- grep (ptn, nms, ignore.case = TRUE)
    }

    if (length (i) != 1) {
        # try initial or terminal punct characters before/after x/y/lon/lat:
        ptn <- paste0 (
            "^[[:punct:]]+", value, "|", value, "[[:punct:]]+$|",
            "^[[:punct:]]+", lonlat, "|", lonlat, "[[:punct:]]+$"
        )
        i <- grep (ptn, nms, ignore.case = TRUE)
    }

    return (i)
}

#' convert_to_matrix
#'
#' Convert a rectangular object to a matrix of two columns, lon and lat.
#'
#' @param obj Rectagular object
#' @return Numeric matrix of two columns: lon and lat
#' @noRd
convert_to_matrix <- function (obj) {

    xy_cols <- find_xy_cols (obj)
    if (is.vector (obj)) {
        obj <- matrix (obj, nrow = 1)
    }
    if (is.numeric (obj)) {

        cbind (obj [, xy_cols [1]], obj [, xy_cols [2]])
    } else {

        if (!(is.numeric (obj [, xy_cols [1], drop = TRUE]) &
            is.numeric (obj [, xy_cols [2], drop = TRUE]))) {

            cbind (
                as.numeric (obj [, xy_cols [1], drop = TRUE]),
                as.numeric (obj [, xy_cols [2], drop = TRUE])
            )
        } else {

            cbind (obj [[xy_cols [1]]], obj [[xy_cols [2]]]) # nolint
        }
    }
}

check_max_d <- function (d, measure) {

    if (max (d, na.rm = TRUE) > 100000) {
        message (
            "Maximum distance is > 100km. The 'cheap' measure is ",
            "inaccurate over such\nlarge distances, you'd likely ",
            "be better using a different 'measure', \n", 
            "one of 'haversine', 'vincenty', or 'geodesic'. "
        )
    }
}

chk_is_num_len_1 <- function (x, xname) {

    if (!is.numeric (x)) {
        stop (xname, " must be numeric")
    }
    if (length (x) > 1) {
        stop (xname, " must be a single value")
    }
}
