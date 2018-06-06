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
        iy <- which (grepl ("y", nms, ignore.case = TRUE) |
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
    if (is.numeric (obj))
        cbind (as.numeric (obj [, xy_cols [1]]),
               as.numeric (obj [, xy_cols [2]]))
    else
        cbind (as.numeric (obj [[xy_cols [1] ]]),
               as.numeric (obj [[xy_cols [2] ]]))
}
