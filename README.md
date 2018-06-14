<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/hypertidy/geodist.svg)](https://travis-ci.org/hypertidy/geodist)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![codecov](https://codecov.io/gh/hypertidy/geodist/branch/master/graph/badge.svg)](https://codecov.io/gh/hypertidy/geodist)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/geodist)](http://cran.r-project.org/web/packages/geodist)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/geodist)

# geodist

An ultra-lightweight, zero-dependency package for very fast calculation
of geodesic distances.

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
is, some variant of `lon`/`lat`, or `x`/`y`).

``` r
library(geodist)
```

``` r
# current verison
packageVersion("geodist")
#> [1] '0.0.0.9000'
```

Input can be in arbtirary rectangular format.

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
x <- cbind (-180 + 360 * runif (n),
             -90 + 100 * runif (n),
             seq (n), runif (n))
colnames (x) <- c ("lon", "lat", "a", "b")
dim (geodist (x))
#> [1] 10 10
```

Distances currently implemented are Haversine, Vincenty (spherical), and
the very fast [mapbox cheap
ruler](https://github.com/mapbox/cheap-ruler-cpp/blob/master/include/mapbox/cheap_ruler.hpp)
(see their [blog
post](https://blog.mapbox.com/fast-geodesic-approximations-with-cheap-ruler-106f229ad016)).
These are determined by the `measure` parameter which takes one of
`"haversine"`, `"vincenty"`, or `"cheap"`.

It is also possible to calculate sequential distances along the rows of
a rentangular object:

``` r
d <- geodist (x, measure = "vincenty", sequential = TRUE)
nrow (x); length (d)
#> [1] 10
#> [1] 9
d <- geodist (x, measure = "vincenty", sequential = TRUE, pad = TRUE)
nrow (x); length (d); head (d)
#> [1] 10
#> [1] 10
#>           [,1]
#> [1,]        NA
#> [2,] 10728.064
#> [3,]  5371.735
#> [4,]  8515.388
#> [5,]  5857.710
#> [6,]  8645.188
```

The `pad` argument pre-pends an `NA` to return a vector commensurate
with input dimensions, and so able to be used directly in piped `dplyr`
operations.

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
#> 2 geodist(x, measure = "vincenty")           10   0.675     1.00
#> 1                     sf_dist(xsf)           10   7.081    10.49
```

The [`geosphere` package](https://cran.r-project.org/package=geosphere)
also offers sequential calculation which is benchmarked with the
following
code:

``` r
fgeodist <- function () geodist (x, measure = "vincenty", sequential = TRUE)
fgeosphere <- function () geosphere::distVincentySphere (x)
rbenchmark::benchmark (replications = 10, order = "test",
                       fgeodist (),
                       fgeosphere ()) [, 1:4]
#>           test replications elapsed relative
#> 1   fgeodist()           10   0.004     1.00
#> 2 fgeosphere()           10   0.053    13.25
```

The [mapbox cheap ruler
algorithm](https://github.com/mapbox/cheap-ruler-cpp) provides
approximate yet very fast distance calculations within small (typically
intra-urban) areas. Speed advantages are demonstrated in the following
code:

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
#> 1     d1 <- geodist(x, measure = "cheap")           10   0.248    1.000
#> 2 d2 <- geodist(x, measure = "haversine")           10   0.301    1.214
#> 3  d3 <- geodist(x, measure = "vincenty")           10   0.393    1.585
```

### Test Results

``` r
require (devtools)
require (testthat)
```

``` r
date()
#> [1] "Thu Jun 14 12:04:52 2018"
devtools::test("tests/")
#> Loading geodist
#> Testing geodist
#> ✔ | OK F W S | Context
✔ | 26       | geodist [0.1 s]
#> 
#> ══ Results ════════════════════════════════════════════════════════════════
#> Duration: 0.2 s
#> 
#> OK:       26
#> Failed:   0
#> Warnings: 0
#> Skipped:  0
```
