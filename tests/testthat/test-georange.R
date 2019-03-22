context("georange")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

test_that ("sequential structure", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              expect_message (d1 <- georange (x, sequential = TRUE,
                                              measure = "haversine"),
                              "object has no named columns")
              colnames (x) <- colnames (y) <- c ("x", "y")
              expect_silent (d1 <- georange (x, sequential = TRUE,
                                             measure = "haversine"))
              expect_message (d2 <- georange (x, y, sequential = TRUE,
                                              measure = "haversine"),
                              "Sequential distances calculated along values of")
              expect_identical (d1, d2)

              expect_equal (length (d1), 2)
              expect_equal (names (d1), c ("minimum", "maximum"))
              expect_true (d1 [2] > d1 [1])
              d3 <- georange (x, sequential = TRUE, measure = "haversine")
              expect_equal (length (d3), 2)
})

test_that ("different measures", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              colnames (x) <- c ("x", "y")
              d1 <- georange (x, sequential = TRUE, measure = "cheap")
              d2 <- georange (x, sequential = TRUE, measure = "haversine")
              d3 <- georange (x, sequential = TRUE, measure = "vincenty")
              d4 <- georange (x, sequential = TRUE, measure = "geodesic")
              expect_true (!identical (d1, d2))
              if (test_all) # haversine and vincenty are sometimes identical
              {
                  expect_true (!identical (d1, d3))
                  expect_true (!identical (d1, d4))
                  expect_true (!identical (d2, d3))
                  expect_true (!identical (d2, d4))
                  expect_true (!identical (d3, d4))
              }
              d5 <- georange (x, sequential = TRUE)
              expect_identical (d1, d5)
              expect_error (d6 <- georange (x, sequential = TRUE,
                                            measure = "junk"),
                            "'arg' should be one of")
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

test_that("range structure for x only", {
              n <- 100
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              colnames (x) <- c ("x", "y")
              d1 <- georange (x, measure = "haversine")
              d2 <- havdist (x)
              diag (d2) <- NA
              d2 <- range (d2, na.rm = TRUE)
              names (d2) <- c ("minimum", "maximum")
              if (test_all)
                  expect_identical (d1, d2)
})

test_that("matrix structure for x y", {
              n <- 100
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- georange (x, y, measure = "haversine")
              d2 <- havdist (x, y)
              d2 <- range (d2, na.rm = TRUE)
              names (d2) <- c ("minimum", "maximum")
              if (test_all)
                  expect_identical (d1, d2)
})

test_that("geodesic extreme cases", {
              # TODO: Write some real tests here
              x <- rbind (c (0, 0),
                          c (0, 1))
              colnames (x) <- c ("x", "y")
              d <- georange (x, measure = "geodesic")
              expect_equal (length (d), 2)

              x <- rbind (c (0, 0),
                          c (0, 90))
              colnames (x) <- c ("x", "y")
              d <- georange (x, measure = "geodesic")
              expect_equal (length (d), 2)

              x <- rbind (c (0, 0),
                          c (1, 0))
              colnames (x) <- c ("x", "y")
              d <- georange (x, measure = "geodesic")
              expect_equal (length (d), 2)

              x <- rbind (c (0, 0),
                          c (180, 0))
              colnames (x) <- c ("x", "y")
              d <- georange (x, measure = "geodesic")
              expect_equal (length (d), 2)
})

test_that ("range measures for x", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              colnames (x) <- c ("x", "y")
              d1 <- georange (x, sequential = TRUE, measure = "haversine")
              d2 <- georange (x, measure = "haversine")
              d3 <- georange (x, measure = "vincenty")
              d4 <- georange (x, measure = "cheap")
              d5 <- georange (x, measure = "geodesic")
              expect_true (!identical (d1, d2))
              expect_true (!identical (d1, d3))
              expect_true (!identical (d1, d4))
              expect_true (!identical (d1, d5))
              expect_true (!identical (d2, d3))
              expect_true (!identical (d2, d4))
              expect_true (!identical (d2, d5))
              expect_true (!identical (d3, d4))
              expect_true (!identical (d3, d5))
              expect_true (!identical (d4, d5))
})

test_that ("range measures for xy", {
              n <- 1e2
              x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d1 <- georange (x, y, measure = "haversine")
              d2 <- georange (x, y, measure = "vincenty")
              d3 <- georange (x, y, measure = "cheap")
              d4 <- georange (x, y, measure = "geodesic")
              expect_true (!identical (d1, d2))
              expect_true (!identical (d1, d3))
              expect_true (!identical (d1, d4))
              expect_true (!identical (d2, d3))
              expect_true (!identical (d2, d4))
              expect_true (!identical (d3, d4))
})
