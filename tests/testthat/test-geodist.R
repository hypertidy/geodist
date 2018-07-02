context("misc tests")

test_all <- identical (Sys.getenv ("MPADGE_LOCAL"), "true")

test_that ("sequential structure", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- geodist (x, sequential = TRUE, measure = "haversine")
              expect_equal (length (d1), nrow (x) - 1)
              d2 <- geodist (x, sequential = TRUE, pad = TRUE,
                             measure = "haversine")
              expect_equal (length (d2), nrow (x))
              # Sequential should equal off-diagonal off full matrix (but note
              # that this test  will fail for "cheap" distances)
              dmat <- geodist (x, measure = "haversine")
              indx <- row (dmat) - col (dmat)
              dmat1 <- split (dmat, indx)["1"][[1]] # first off-diagonal
              expect_identical (d1, dmat1)
              expect_message (d3 <- geodist (x, y, sequential = TRUE),
                              "Sequential distances calculated along values")
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
              d1 <- geodist (x, measure = "haversine")
              d2 <- havdist (x)
              if (test_all)
                  expect_identical (d1, d2)
})

test_that("matrix structure for x y", {
              n <- 100
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- geodist (x, y, measure = "haversine")
              d2 <- havdist (x, y)
              if (test_all)
                  expect_identical (d1, d2)
})

test_that("geodesic extreme cases", {
              # TODO: Write some real tests here
              x <- rbind (c (0, 0),
                          c (0, 1))
              colnames (x) <- c ("x", "y")
              d <- geodist (x, measure = "geodesic")
              expect_true (sum (diag (d)) == 0)

              x <- rbind (c (0, 0),
                          c (0, 90))
              colnames (x) <- c ("x", "y")
              d <- geodist (x, measure = "geodesic")
              expect_true (sum (diag (d)) == 0)

              x <- rbind (c (0, 0),
                          c (1, 0))
              colnames (x) <- c ("x", "y")
              d <- geodist (x, measure = "geodesic")
              expect_true (sum (diag (d)) == 0)

              x <- rbind (c (0, 0),
                          c (180, 0))
              colnames (x) <- c ("x", "y")
              d <- geodist (x, measure = "geodesic")
              expect_true (sum (diag (d)) == 0)
              m <- 20003930
              expect_true (abs (d [1, 2] - m) < 2) # it doesn't equal zero
              expect_true (abs (d [2, 1] - m) < 2) # it doesn't equal zero
})

test_that ("geodist_benchmark", {
               d <- geodist_benchmark (lat = 1, d = 100, n = 100)
               expect_is (d, "matrix")
               expect_equal (nrow (d), 2)
               expect_equal (ncol (d), 3)
               expect_equal (rownames (d), c ("absolute", "relative"))
               expect_equal (colnames (d),
                             c ("haversine", "vincenty", "cheap"))
})
