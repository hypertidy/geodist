context("geodist measures")

# all measures compared against geodesic

test_that("measures", {
              n <- 1e2
              #x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
              #y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
              x <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
              y <- cbind (-10 + 20 * runif (2 * n), -10 + 20 * runif (2 * n))
              colnames (x) <- colnames (y) <- c ("x", "y")
              d0_x <- geodist (x, measure = "geodesic")
              d0_xy <- geodist (x, y, measure = "geodesic")
              d0_seq <- geodist (x, measure = "geodesic", sequential = TRUE)

              measures <- c ("haversine", "vincenty", "vincenty_ellips",
                             "cheap")
              for (m in measures)
              {
                  d1 <- geodist (x, measure = m)
                  expect_true (max (abs (d1 - d0_x)) > 0)
                  expect_true (cor (as.vector (d0_x), as.vector (d1)) > 0.99)
                  d1 <- geodist (x, y, measure = m)
                  expect_true (max (abs (d1 - d0_xy)) > 0)
                  expect_true (cor (as.vector (d0_xy), as.vector (d1)) > 0.99)
                  d1 <- geodist (x, measure = m, sequential = TRUE)
                  expect_true (max (abs (d1 - d0_seq)) > 0)
                  expect_true (cor (as.vector (d0_seq), as.vector (d1)) > 0.99)
              }
})
