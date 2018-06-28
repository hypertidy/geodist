context("geodist")
# tibble not tested to avoid pacakge suggests, but just needs these lines
#x <- tibble::tibble (x = -180 + 360 * runif (n),
#                     y = -90 + 180 * runif (n))
#y <- tibble::tibble (x = -180 + 360 * runif (2 * n),
#                     y = -90 + 180 * runif (2 * n))

test_that("geodist with df", {
              n <- 1e2
              x <- data.frame (x = -180 + 360 * runif (n),
                               y = -90 + 180 * runif (n))
              y <- data.frame (x = -180 + 360 * runif (2 * n),
                               y = -90 + 180 * runif (2 * n))
              d1 <- geodist (x)
              expect_equal (dim (d1), c (n, n))
              diag (d1) <- Inf
              expect_true (all (d1 >= 0))

              d2 <- geodist (x, y)
              expect_equal (dim (d2), c (n, 2 * n))
              diag (d2) <- Inf
              expect_true (all (d2 >= 0))
})

test_that("geodist with matrix", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              expect_message (d1 <- geodist (x), "object has no named columns")
              colnames (x) <- c ("x", "y")
              expect_silent (d1 <- geodist (x))
              expect_equal (dim (d1), c (n, n))
              diag (d1) <- Inf
              expect_true (all (d1 >= 0))

              expect_message (d2 <- geodist (x, y), "object has no named columns")
              expect_equal (dim (d2), c (n, 2 * n))
              diag (d2) <- Inf
              expect_true (all (d2 > 0))
              expect_true (all (d2 >= 0))
})

test_that("haversine and vincenty", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- geodist (x, measure = "haversine")
              d2 <- geodist (x, measure = "vincenty")
              expect_true (!identical (d1, d2))
              d1 <- geodist (x, y, measure = "haversine")
              d2 <- geodist (x, y, measure = "vincenty")
              expect_true (!identical (d1, d2))
})

test_that("haversine and cheap", {
              n <- 1e2
              dx <- dy <- 0.01
              x <- cbind (-100 + dx * runif (n), 20 + dy * runif (n))
              y <- cbind (-100 + dx * runif (2 * n), 20 + dy * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- geodist (x, measure = "haversine")
              d2 <- geodist (x, measure = "cheap")
              expect_true (!identical (d1, d2))
              d1 <- geodist (x, y, measure = "haversine")
              d2 <- geodist (x, y, measure = "cheap")
              expect_true (!identical (d1, d2))
})

test_that("haversine and geodesic", {
              n <- 1e2
              dx <- dy <- 0.01
              x <- cbind (-100 + dx * runif (n), 20 + dy * runif (n))
              y <- cbind (-100 + dx * runif (2 * n), 20 + dy * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- geodist (x, measure = "haversine")
              d2 <- geodist (x, measure = "geodesic")
              expect_true (!identical (d1, d2))
              d1 <- geodist (x, y, measure = "haversine")
              d2 <- geodist (x, y, measure = "geodesic")
              expect_true (!identical (d1, d2))
})

test_that("vincenty and geodesic", {
              n <- 1e2
              dx <- dy <- 0.01
              x <- cbind (-100 + dx * runif (n), 20 + dy * runif (n))
              y <- cbind (-100 + dx * runif (2 * n), 20 + dy * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- geodist (x, measure = "vincenty")
              d2 <- geodist (x, measure = "geodesic")
              expect_true (!identical (d1, d2))
              d1 <- geodist (x, y, measure = "vincenty")
              d2 <- geodist (x, y, measure = "geodesic")
              expect_true (!identical (d1, d2))
})

test_that ("sequential structure", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- geodist (x, sequential = TRUE)
              expect_equal (length (d1), nrow (x) - 1)
              d2 <- geodist (x, sequential = TRUE, pad = TRUE)
              expect_equal (length (d2), nrow (x))
              # Sequential should equal off-diagonal off full matrix:
              dmat <- geodist (x)
              indx <- row (dmat) - col (dmat)
              dmat1 <- split (dmat, indx)["1"][[1]] # first off-diagonal
              expect_identical (d1, dmat1)
              expect_message (d3 <- geodist (x, y, sequential = TRUE),
                              "Sequential distances calculated along values")
})

test_that ("sequential measures", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              colnames (x) <- c ("x", "y")
              d1 <- geodist (x, sequential = TRUE, measure = "haversine")
              d2 <- geodist (x, sequential = TRUE, measure = "vincenty")
              d3 <- geodist (x, sequential = TRUE, measure = "cheap")
              d4 <- geodist (x, sequential = TRUE, measure = "geodesic")
              expect_error (d4 <- geodist (x, sequential = TRUE,
                                           measure = "blah"))
              expect_true (max (abs (d1 - d2)) > 0)
              expect_true (max (abs (d1 - d3)) > 0)
              expect_true (max (abs (d1 - d4)) > 0)
              expect_true (max (abs (d2 - d3)) > 0)
              expect_true (max (abs (d2 - d4)) > 0)
              expect_true (max (abs (d3 - d4)) > 0)
})

havdist <- function (x, y)
{
    if (missing (y))
        y <- x
    x1mat <- array (x [, 1], dim = c (nrow (x), nrow (y)))
    x2mat <- array (x [, 2], dim = c (nrow (x), nrow (y)))
    y1mat <- t (array (y [, 1], dim = c (nrow (y), nrow (x))))
    y2mat <- t (array (y [, 2], dim = c (nrow (y), nrow (x))))
    xd <- (x1mat - y1mat) * pi / 180
    yd <- (x2mat - y2mat) * pi / 180
    sxd <- sin (xd / 2)
    syd <- sin (yd / 2)
    earth <- 6378137 # radius of Earth in m
    d <- syd * syd +
        cos (x2mat * pi / 180) * cos (y2mat * pi / 180) * sxd * sxd
    2 * earth * asin (sqrt (d));
}

test_that("matrix structure for x only", {
              n <- 100
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              colnames (x) <- c ("x", "y")
              d1 <- geodist (x)
              d2 <- havdist (x)
              expect_identical (d1, d2)
})

test_that("matrix structure for x y", {
              n <- 100
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- geodist (x, y)
              d2 <- havdist (x, y)
              expect_identical (d1, d2)
})

test_that("geodesic extreme cases", {
              x <- rbind (c (0, 0),
                          c (0, 1))
              colnames (x) <- c ("x", "y")
              d <- geodist (x, measure = "geodesic")
              expect_true (sum (diag (d)) > 0)

              x <- rbind (c (0, 0),
                          c (0, 90))
              colnames (x) <- c ("x", "y")
              d <- geodist (x, measure = "geodesic")
              expect_true (sum (diag (d)) > 0)

              x <- rbind (c (0, 0),
                          c (1, 0))
              colnames (x) <- c ("x", "y")
              d <- geodist (x, measure = "geodesic")
              expect_true (sum (diag (d)) > 0)

              x <- rbind (c (0, 0),
                          c (180, 0))
              colnames (x) <- c ("x", "y")
              d <- geodist (x, measure = "geodesic")
              expect_true (sum (diag (d)) > 0)
              m <- 20003930
              expect_true (abs (d [1, 2] - m) < 2) # it doesn't equal zero
              expect_true (abs (d [2, 1] - m) < 2) # it doesn't equal zero
})
