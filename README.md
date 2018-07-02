<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/hypertidy/geodist.svg)](https://travis-ci.org/hypertidy/geodist)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/hypertidy/geodist?branch=master&svg=true)](https://ci.appveyor.com/project/hypertidy/geodist)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/hypertidy/geodist/branch/master/graph/badge.svg)](https://codecov.io/gh/hypertidy/geodist)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/geodist)](http://cran.r-project.org/web/packages/geodist)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/geodist)

# geodist

An ultra-lightweight, zero-dependency package for very fast calculation
of geodesic distances. Main eponymous function, `geodist()`, accepts
only one or two primary arguments, which must be rectagular objects with
unambiguously labelled longitude and latitude columns (that is, some
variant of `lon`/`lat`, or `x`/`y`).

``` r
n <- 50
x <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
y <- cbind (-10 + 20 * runif (2 * n), -10 + 20 * runif (2 * n))
colnames (x) <- colnames (y) <- c ("x", "y")
d0 <- geodist (x) # A 50-by-50 matrix
d1 <- geodist (x, y) # A 50-by-100 matrix
d2 <- geodist (x, sequential = TRUE) # Vector of length 49
d2 <- geodist (x, sequential = TRUE, pad = TRUE) # Vector of length 50
```

## Installation

You can install geodist from github with:

``` r
# install.packages("devtools")
devtools::install_github("hypertidy/geodist")
```

``` r
library(geodist)
```

``` r
# current verison
packageVersion("geodist")
#> [1] '0.0.1'
```

## Detailed Usage

Input(s) to the `geodist()` function can be in arbitrary rectangular
format.

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

Distances currently implemented are Haversine, Vincenty (spherical and
elliptical)), the very fast [mapbox cheap
ruler](https://github.com/mapbox/cheap-ruler-cpp/blob/master/include/mapbox/cheap_ruler.hpp)
(see their [blog
post](https://blog.mapbox.com/fast-geodesic-approximations-with-cheap-ruler-106f229ad016)),
and the “reference” implementation of [Karney
(2013)](https://link.springer.com/content/pdf/10.1007/s00190-012-0578-z.pdf),
as implemented in the package
[`sf`](https://cran.r-project.org/package=sf). (Note that `geodist` does
not accept [`sf`](https://cran.r-project.org/package=sf)-format objects;
the [`sf`](https://cran.r-project.org/package=sf) package itself should
be used for that.) The [mapbox cheap ruler
algorithm](https://github.com/mapbox/cheap-ruler-cpp) is intended to
provide approximate yet very fast distance calculations within small
areas (tens to a few hundred kilometres across).

### Benchmarks of geodetic accuracy

The `geodist_benchmark()` function - the only other function provided by
the `geodist` package - compares the accuracy of the different metrics
to the nanometre-accuracy standard of [Karney
(2013)](https://link.springer.com/content/pdf/10.1007/s00190-012-0578-z.pdf).

``` r
geodist_benchmark (lat = 30, d = 1000)
#>            haversine    vincenty       cheap
#> absolute 0.821979685 0.821979685 0.573589772
#> relative 0.002206036 0.002206036 0.001613667
```

All distances (`d)` are in metres, and all measures are accurate to
within 1m over distances out to several km (at the chosen latitude of 30
degrees). The following plots compare the absolute and relative
accuracies of the different distance measures implemented here. The
mapbox cheap ruler algorithm is the most accurate for distances out to
around 100km, beyond which it becomes extremely inaccurate. Average
relative errors of Vincenty distances remain generally constant at
around 0.2%, while relative errors of cheap-ruler distances out to 100km
are around 0.16%.

![](vignettes/fig1.png)

### Performance comparison

The following code demonstrates the relative speed advantages of the
different distance measures implemented in the `geodist` package.

``` r
n <- 1e3
dx <- dy <- 0.01
x <- cbind (-100 + dx * runif (n), 20 + dy * runif (n))
y <- cbind (-100 + dx * runif (2 * n), 20 + dy * runif (2 * n))
colnames (x) <- colnames (y) <- c ("x", "y")
rbenchmark::benchmark (replications = 10, order = "test",
                       d1 <- geodist (x, measure = "cheap"),
                       d2 <- geodist (x, measure = "haversine"),
                       d3 <- geodist (x, measure = "vincenty"),
                       d4 <- geodist (x, measure = "geodesic")) [, 1:4]
#>                                      test replications elapsed relative
#> 1     d1 <- geodist(x, measure = "cheap")           10   0.158    1.000
#> 2 d2 <- geodist(x, measure = "haversine")           10   0.243    1.538
#> 3  d3 <- geodist(x, measure = "vincenty")           10   0.394    2.494
#> 4  d4 <- geodist(x, measure = "geodesic")           10   4.559   28.854
```

Geodesic distance calculation is available in the [`sf`
package](https://cran.r-project.org/package=sf). Comparing computation
speeds requires conversion of sets of numeric lon-lat points to `sf`
form with the following code:

``` r
require (magrittr)
#> Loading required package: magrittr
x_to_sf <- function (x)
{
    sapply (seq (nrow (x)), function (i)
            sf::st_point (x [i, ]) %>%
                sf::st_sfc ()) %>%
    sf::st_sfc (crs = 4326)
}
```

``` r
n <- 1e2
x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
colnames (x) <- c ("x", "y")
xsf <- x_to_sf (x)
sf_dist <- function (xsf) sf::st_distance (xsf, xsf)
geo_dist <- function (x) geodist (x, measure = "geodesic")
rbenchmark::benchmark (replications = 10, order = "test",
                      sf_dist (xsf),
                      geo_dist (x)) [, 1:4]
#> Linking to GEOS 3.5.0, GDAL 2.1.3, proj.4 4.9.2
#>           test replications elapsed relative
#> 2  geo_dist(x)           10   0.080    1.000
#> 1 sf_dist(xsf)           10   0.353    4.412
```

Confirm that the two give almost identical results:

``` r
ds <- matrix (as.numeric (sf_dist (xsf)), nrow = length (xsf))
dg <- geodist (x, measure = "geodesic")
formatC (max (abs (ds - dg)), format = "e")
#> [1] "7.4506e-09"
```

All results are in metres, so the two differ by only around 10
nanometres.

The [`geosphere` package](https://cran.r-project.org/package=geosphere)
also offers sequential calculation which is benchmarked with the
following
code:

``` r
fgeodist <- function () geodist (x, measure = "vincenty", sequential = TRUE)
fgeosph <- function () geosphere::distVincentySphere (x)
rbenchmark::benchmark (replications = 10, order = "test",
                       fgeodist (),
                       fgeosph ()) [, 1:4]
#>         test replications elapsed relative
#> 1 fgeodist()           10   0.024    1.000
#> 2  fgeosph()           10   0.100    4.167
```

`geodist` is thus around 4 times faster than both `sf` for highly accurate
geodesic distance calculations, and `geosphere` for calculation of sequential
distances.

### Test Results

``` r
require (devtools)
require (testthat)
```

``` r
date()
#> [1] "Mon Jul  2 11:21:57 2018"
devtools::test("tests/")
#> Loading geodist
#> Testing geodist
#> ✔ | OK F W S | Context
✔ | 17       | misc tests [0.1 s]
✔ | 12       | geodist input formats
✔ | 18       | geodist measures
#> 
#> ══ Results ════════════════════════════════════════════════════════════════
#> Duration: 0.3 s
#> 
#> OK:       47
#> Failed:   0
#> Warnings: 0
#> Skipped:  0
```
