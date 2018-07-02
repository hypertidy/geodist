context("geodist measures")

test_all <- identical (Sys.getenv ("MPADGE_LOCAL"), "true")

# all measures compared against geodesic
cor_lim <- 0.9

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

              measures <- c ("haversine", "vincenty", "cheap")
              for (m in measures)
              {
                  d1_x <- geodist (x, measure = m)
                  d1_xy <- geodist (x, y, measure = m)
                  d1_seq <- geodist (x, measure = m, sequential = TRUE)

                  expect_true (max (abs (d1_x - d0_x)) > 0)
                  expect_true (max (abs (d1_xy - d0_xy)) > 0)
                  expect_true (max (abs (d1_seq - d0_seq)) > 0)

                  if (test_all)
                  {
                      expect_true (cor (as.vector (d0_x),
                                        as.vector (d1_x)) > cor_lim)
                      expect_true (cor (as.vector (d0_xy),
                                        as.vector (d1_xy)) > cor_lim)
                      expect_true (cor (as.vector (d0_seq),
                                        as.vector (d1_seq)) > cor_lim)
                  }
              }
})
