test_all <- identical (Sys.getenv ("MPADGE_LOCAL"), "true")

test_that ("geodist min", {

    n <- 1e2
    x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
    y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
    colnames (x) <- colnames (y) <- c ("x", "y")

    d0 <- geodist (x, y, measure = "Haversine")
    index0 <- apply (d0, 1, which.min)
    index1 <- geodist_min (x, y, measure = "Haversine")
    expect_identical (index0, index1)

    d0 <- geodist (x, y, measure = "cheap")
    index0 <- apply (d0, 1, which.min)
    index1 <- geodist_min (x, y, measure = "cheap")
    expect_identical (index0, index1)

    d0 <- geodist (x, y, measure = "vincenty")
    index0 <- apply (d0, 1, which.min)
    index1 <- geodist_min (x, y, measure = "vincenty")
    expect_identical (index0, index1)

    d0 <- geodist (x, y, measure = "geodesic")
    index0 <- apply (d0, 1, which.min)
    index1 <- geodist_min (x, y, measure = "geodesic")
    expect_identical (index0, index1)
})
