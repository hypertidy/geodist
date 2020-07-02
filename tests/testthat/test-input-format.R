context("geodist input formats")
# tibble not tested to avoid pacakge suggests, but just needs these lines
#x <- tibble::tibble (x = -180 + 360 * runif (n),
#                     y = -90 + 180 * runif (n))
#y <- tibble::tibble (x = -180 + 360 * runif (2 * n),
#                     y = -90 + 180 * runif (2 * n))

test_all <- identical (Sys.getenv ("MPADGE_LOCAL"), "true")

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
              if (test_all)
                  expect_message (d1 <- geodist (x),
                                  "Maximum distance is > 100km")
              else
                  d1 <- geodist (x)
              expect_equal (dim (d1), c (n, n))
              diag (d1) <- Inf
              expect_true (all (d1 >= 0))

              expect_message (d2 <- geodist (x, y), "object has no named columns")
              expect_equal (dim (d2), c (n, 2 * n))
              diag (d2) <- Inf
              expect_true (all (d2 > 0))
              expect_true (all (d2 >= 0))
})

test_that ("other columns", {
               n <- 50
               x <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
               y <- cbind (-10 + 20 * runif (2 * n), -10 + 20 * runif (2 * n))
               colnames (x) <- colnames (y) <- c ("x", "y")
               d <- geodist (x, y)
               x1 <- cbind (newx = runif (n), newy = runif (n), x)
               y1 <- cbind (newx = runif (n), newy = runif (n), y)
               d1 <- geodist (x1, y1)
               expect_identical (d, d1)
})

test_that ("column names, ", {
               n <- 50
               x <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
               y <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
               colnames (x) <- colnames (y) <- c ("x", "x")
               expect_error (d <- geodist (x, y, paired = TRUE),
                             paste0 ("Unable to determine longitude and ",
                                     "latitude columns"))
               colnames (x) <- colnames (y) <- c ("_x_x_", "_x_y_")
               expect_error (d <- geodist (x, y, paired = TRUE),
                             paste0 ("Unable to determine longitude and ",
                                     "latitude columns"))
               colnames (x) <- colnames (y) <- c ("_x_x", "_x_y")
               if (test_all)
                   expect_message (d1 <- geodist (x, y, paired = TRUE),
                                  "Maximum distance is > 100km")
               else
                   d1 <- geodist (x, y, paired = TRUE)

               colnames (x) <- colnames (y) <- c ("_a_x_", "_a_y_")
               if (test_all)
                   expect_message (d2 <- geodist (x, y, paired = TRUE),
                                   "Maximum distance is > 100km")
               else
                   d2 <- geodist (x, y, paired = TRUE)
               expect_identical (d1, d2)
})
