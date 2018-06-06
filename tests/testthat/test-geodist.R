context("geodist")
test_that("geodist with tibble", {
              n <- 1e2
              x <- tibble::tibble (x = -180 + 360 * runif (n),
                                   y = -90 + 180 * runif (n))
              y <- tibble::tibble (x = -180 + 360 * runif (2 * n),
                                   y = -90 + 180 * runif (2 * n))
              d1 <- geodist (x)
              expect_equal (dim (d1), c (n, n))
              expect_true (all (d1 >= 0))

              d2 <- geodist (x, y)
              expect_equal (dim (d2), c (n, 2 * n))
              expect_true (all (d2 >= 0))
})

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
