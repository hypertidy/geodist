<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/hypertidy/geodist.svg)](https://travis-ci.org/hypertidy/geodist)
[![Project Status: Concept - Minimal or no implementation has been done
yet.](http://www.repostatus.org/badges/0.1.0/concept.svg)](http://www.repostatus.org/#concept)
[![codecov](https://codecov.io/gh/hypertidy/geodist/branch/master/graph/badge.svg)](https://codecov.io/gh/hypertidy/geodist)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/geodist)](http://cran.r-project.org/web/packages/geodist)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/geodist)

# geodist

An ultra-lightweight, zero-dependency package for calculation of
geodesic distances.

## Installation

You can install geodist from github with:

``` r
# install.packages("devtools")
devtools::install_github("hypertidy/geodist")
```

## Usage

`geodist` contains only one eponymous function which accepts only one or
two primary arguments, each of which must be some kind of rectangular
object with unambiguously labelled longitude and latitude columns (that
is, some variant of `lon`/`lat`, or `x`/`y`). The only other argument is
`measure`, specifying one of Haversine or Vincenty great-circle
distances.

    #> Loading geodist

``` r
library(geodist)
```

``` r
# current verison
packageVersion("geodist")
#> [1] '0.0.0.9000'
```

``` r
n <- 1e1
x <- tibble::tibble (x = -180 + 360 * runif (n),
                     y = -90 + 180 * runif (n))
dim (geodist (x))
#> [1] 10 10
y <- tibble::tibble (x = -180 + 360 * runif (2 * n),
                     y = -90 + 180 * runif (2 * n))
dim (geodist (x, y))
#> [1] 10 20
```

### Comparison

Geodesic distance calculation is available in the [`sf`
package](https://cran.r-project.org/package=sf), requiring only
conversion of sets of numeric lon-lat points to `sf` form with the
following code:

``` r
require (magrittr)
#> Loading required package: magrittr
x_to_sf <- function (x)
{
    sapply (seq (nrow (x)), function (i)
            sf::st_point (x [i, ]) %>%
                sf::st_sfc ()) %>%
    sf::st_sfc ()
}
```

Note that the relevant function, `sf::st_distance()` calculates the
Vicenty distance via `liblwgeom`, so `measure = "vincenty"` is used in
the following comparison.

``` r
n <- 1e3
x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
colnames (x) <- c ("x", "y")
xsf <- x_to_sf (x)
sf_dist <- function (x) sf::st_distance (x, x)
rbenchmark::benchmark (replications = 10, order = "test",
                      sf_dist (xsf),
                      geodist (x, measure = "vincenty")) [, 1:4]
#>                               test replications elapsed relative
#> 2 geodist(x, measure = "vincenty")           10   0.664    1.000
#> 1                     sf_dist(xsf)           10   6.990   10.527
```

The `geodist` package also implements the [mapbox cheap ruler
algorithm](https://github.com/mapbox/cheap-ruler-cpp) (described in this
[blog
post](https://blog.mapbox.com/fast-geodesic-approximations-with-cheap-ruler-106f229ad016),
providing approximate yet very fast distance calculations within small
(typically intra-urban) areas. Speed advantages are demonstrated in the
following code:

``` r
n <- 1e3
dx <- dy <- 0.01
x <- cbind (-100 + dx * runif (n), 20 + dy * runif (n))
y <- cbind (-100 + dx * runif (2 * n), 20 + dy * runif (2 * n))
colnames (x) <- colnames (y) <- c ("x", "y")
rbenchmark::benchmark (replications = 10, order = "test",
                       d1 <- geodist (x, measure = "cheap"),
                       d2 <- geodist (x, measure = "haversine"),
                       d3 <- geodist (x, measure = "vincenty")) [, 1:4]
#>                                      test replications elapsed relative
#> 1     d1 <- geodist(x, measure = "cheap")           10   0.244    1.000
#> 2 d2 <- geodist(x, measure = "haversine")           10   0.279    1.143
#> 3  d3 <- geodist(x, measure = "vincenty")           10   0.379    1.553
```

### Test Results

``` r
require (devtools)
require (testthat)
```

``` r
date()
#> [1] "Wed Jun 13 17:28:59 2018"
devtools::test("tests/")
#> Loading geodist
#> Testing geodist
#> ✔ | OK F W S | Context
#> 
⠏ |  0       | geodist
⠋ |  1       | geodist
⠙ |  2       | geodist
⠹ |  3       | geodist
⠸ |  4       | geodist
⠼ |  5       | geodist
⠴ |  6       | geodist
⠦ |  7       | geodist
⠧ |  8       | geodist
⠇ |  9       | geodist
⠏ | 10       | geodist
⠋ | 11       | geodist
⠙ | 12       | geodist
⠹ | 13       | geodist
⠸ | 14       | geodist
⠼ | 15       | geodist
⠴ | 16       | geodist
⠦ | 17       | geodist
⠧ | 18       | geodist
✔ | 18       | geodist [0.1 s]
#> 
#> ══ Results ════════════════════════════════════════════════════════════════
#> Duration: 0.2 s
#> 
#> OK:       18
#> Failed:   0
#> Warnings: 0
#> Skipped:  0
#> 
#> Way to go!
```
