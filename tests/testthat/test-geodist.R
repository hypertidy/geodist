test_all <- identical (Sys.getenv ("MPADGE_LOCAL"), "true")

test_that ("sequential structure", {
    n <- 1e2
    x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
    y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
    colnames (x) <- colnames (y) <- c ("x", "y")
    d1 <- geodist (x, sequential = TRUE, measure = "haversine")
    expect_equal (length (d1), nrow (x) - 1)
    d2 <- geodist (x,
        sequential = TRUE, pad = TRUE,
        measure = "haversine"
    )
    expect_equal (length (d2), nrow (x))
    # Sequential should equal off-diagonal of full matrix (but note
    # that this test  will fail for "cheap" distances)
    dmat <- geodist (x, measure = "haversine")
    indx <- row (dmat) - col (dmat)
    dmat1 <- split (dmat, indx) ["1"] [[1]] # first off-diagonal
    expect_identical (d1, dmat1)
    expect_message (
        d3 <- geodist (x, y, sequential = TRUE),
        "Sequential distances calculated along values"
    )
})

havdist <- function (x, y) {

    if (missing (y)) {
        y <- x
    }
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
    2 * earth * asin (sqrt (d))
}

test_that ("matrix structure for x only", {
    n <- 100
    x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
    colnames (x) <- c ("x", "y")
    d1 <- geodist (x, measure = "haversine")
    d2 <- havdist (x)
    if (test_all) {
        expect_identical (d1, d2)
    }
})

test_that ("matrix structure for x y", {
    n <- 100
    x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
    y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
    colnames (x) <- colnames (y) <- c ("x", "y")
    d1 <- geodist (x, y, measure = "haversine")
    d2 <- havdist (x, y)
    if (test_all) {
        expect_identical (d1, d2)
    }
})

test_that ("geodesic extreme cases", {
    # TODO: Write some real tests here
    x <- rbind (
        c (0, 0),
        c (0, 1)
    )
    colnames (x) <- c ("x", "y")
    d <- geodist (x, measure = "geodesic")
    expect_true (sum (diag (d)) == 0)

    x <- rbind (
        c (0, 0),
        c (0, 90)
    )
    colnames (x) <- c ("x", "y")
    d <- geodist (x, measure = "geodesic")
    expect_true (sum (diag (d)) == 0)

    x <- rbind (
        c (0, 0),
        c (1, 0)
    )
    colnames (x) <- c ("x", "y")
    d <- geodist (x, measure = "geodesic")
    expect_true (sum (diag (d)) == 0)

    x <- rbind (
        c (0, 0),
        c (180, 0)
    )
    colnames (x) <- c ("x", "y")
    d <- geodist (x, measure = "geodesic")
    expect_true (sum (diag (d)) == 0)
    m <- 20003930
    expect_true (abs (d [1, 2] - m) < 2) # it doesn't equal zero
    expect_true (abs (d [2, 1] - m) < 2) # it doesn't equal zero
})

test_that ("geodist_benchmark", {
    d <- geodist_benchmark (lat = 1, d = 100, n = 100)
    expect_true (inherits (d, "matrix"))
    expect_equal (nrow (d), 2)
    expect_equal (ncol (d), 3)
    expect_equal (rownames (d), c ("absolute", "relative"))
    expect_equal (
        colnames (d),
        c ("haversine", "vincenty", "cheap")
    )

    # benchmarking is restricted to 2 <= n <= 1e3
    expect_error (
        geodist_benchmark (n = 1),
        "Comparisons require at least n = 2 objects"
    )
    expect_error (
        geodist_benchmark (n = 1e3 + 1),
        "nothing to be gained by extending beyond a million"
    )
})

test_that ("geodist paired", {
    n <- 1e2
    x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
    y <- cbind (-180 + 360 * runif (2 * n), -90 + 180 * runif (2 * n))
    colnames (x) <- colnames (y) <- c ("x", "y")
    expect_error (
        d1 <- geodist (x, y, paired = TRUE),
        "x and y must have the same number of rows"
    )
    y <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
    colnames (y) <- c ("x", "y")
    if (test_all) {
        expect_message (
            d1 <- geodist (x, y, paired = TRUE),
            "Maximum distance is > 100km"
        )
    } else {
        d1 <- geodist (x, y, paired = TRUE)
    }
    expect_type (d1, "double")
    expect_equal (length (d1), n)

    expect_silent (d2 <- geodist (x, y,
        paired = TRUE,
        measure = "haversine"
    ))
    expect_silent (d3 <- geodist (x, y,
        paired = TRUE,
        measure = "vincenty"
    ))
    expect_silent (d4 <- geodist (x, y,
        paired = TRUE,
        measure = "geodesic"
    ))
    expect_true (cor (d2, d3) > 0.99)
    expect_true (cor (d2, d4) > 0.99)
})

test_that ("geodist_vec", {
    n <- 1e2
    x1 <- runif (n, -0.1, 0.1)
    y1 <- runif (n, -0.1, 0.1)
    x2 <- runif (n, -0.1, 0.1)
    y2 <- runif (n, -0.1, 0.1)
    x <- cbind ("x" = x1, "y" = y1)
    y <- cbind ("x" = x2, "y" = y2)

    # errors:
    expect_message (
        d <- geodist_vec (x1, y1, x2, y2,
            sequential = TRUE
        ),
        paste0 (
            "Sequential distances calculated along ",
            "values of 'x' only"
        )
    )
    expect_error (
        d <- geodist_vec (),
        "x1 and y1 must be provided"
    )
    expect_error (
        d <- geodist_vec (x1 = x, y1 = y),
        "geodist_vec only accepts vector inputs"
    )
    expect_error (
        d <- geodist_vec (x1, y1 [1:10]),
        "x1 and y1 must have the same length"
    )

    measures <- c ("cheap", "haversine", "vincenty", "geodesic")
    for (m in measures) {

        d1 <- geodist (x, y, paired = TRUE, measure = m)
        d2 <- geodist_vec (x1, y1, x2, y2,
            paired = TRUE, measure = m
        )
        expect_identical (d1, d2)

        d1 <- geodist (x, sequential = TRUE, measure = m)
        d2 <- geodist_vec (x1, y1, sequential = TRUE, measure = m)
        expect_identical (d1, d2)

        d1 <- geodist (x, measure = m)
        d2 <- geodist_vec (x1, y1, measure = m)
        expect_identical (d1, d2)

        d1 <- geodist (x, y, measure = m)
        d2 <- geodist_vec (x1, y1, x2, y2, measure = m)
        expect_identical (d1, d2)
    }
})

test_that ("argument quiet for measure cheap works", {
    x1 <- c (0, 1)
    y1 <- c (0, 1)
    x2 <- c (2, 3)
    y2 <- c (2, 3)
    x <- cbind ("x" = x1, "y" = y1)
    y <- cbind ("x" = x2, "y" = y2)

    expect_message (
        d1 <- geodist (x, y, paired = TRUE),
        "Maximum distance is > 100km"
    )
    expect_silent (d1 <- geodist (x, y, paired = TRUE, quiet = TRUE))

    expect_message (
        d1 <- geodist_vec (x1, y1, x2, y2, paired = TRUE),
        "Maximum distance is > 100km"
    )
    expect_silent (d1 <- geodist_vec (x1, y1, x2, y2,
        paired = TRUE, quiet = TRUE
    ))

})
