<!-- README.md is generated from README.Rmd. Please edit that file -->

[![R build
status](https://github.com/hypertidy/geodist/workflows/R-CMD-check/badge.svg)](https://github.com/hypertidy/geodist/actions?query=workflow%3AR-CMD-check)
[![pkgcheck](https://github.com/hypertidy/geodist/workflows/pkgcheck/badge.svg)](https://github.com/hypertidy/geodist/actions?query=workflow%3Apkgcheck)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/hypertidy/geodist/branch/master/graph/badge.svg)](https://codecov.io/gh/hypertidy/geodist)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/geodist)](http://cran.r-project.org/web/packages/geodist)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/geodist)

# geodist

An ultra-lightweight, zero-dependency package for very fast calculation
of geodesic distances. Main eponymous function, `geodist()`, accepts
only one or two primary arguments, which must be rectangular objects
with unambiguously labelled longitude and latitude columns (that is,
some variant of `lon`/`lat`, or `x`/`y`).

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

You can install latest stable version of `geodist` from CRAN with:

``` r
install.packages ("geodist") # current CRAN version
```

Alternatively, current development versions can be installed using any
of the following options:

``` r
# install.packages("remotes")
remotes::install_git ("https://git.sr.ht/~mpadge/geodist")
remotes::install_git ("https://codeberg.org/hypertidy/geodist")
remotes::install_bitbucket ("hypertidy/geodist")
remotes::install_gitlab ("hypertidy/geodist")
remotes::install_github ("hypertidy/geodist")
```

Then load with

``` r
library (geodist)
packageVersion ("geodist")
#> [1] '0.0.8.6'
```

## Detailed Usage

Input(s) to the `geodist()` function can be in arbitrary rectangular
format.

``` r
n <- 1e1
x <- tibble::tibble (
    x = -180 + 360 * runif (n),
    y = -90 + 180 * runif (n)
)
dim (geodist (x))
#> Maximum distance is > 100km. The 'cheap' measure is inaccurate over such
#> large distances, you'd likely be better using a different 'measure'.
#> [1] 10 10
y <- tibble::tibble (
    x = -180 + 360 * runif (2 * n),
    y = -90 + 180 * runif (2 * n)
)
dim (geodist (x, y))
#> Maximum distance is > 100km. The 'cheap' measure is inaccurate over such
#> large distances, you'd likely be better using a different 'measure'.
#> [1] 10 20
x <- cbind (
    -180 + 360 * runif (n),
    -90 + 100 * runif (n),
    seq (n), runif (n)
)
colnames (x) <- c ("lon", "lat", "a", "b")
dim (geodist (x))
#> Maximum distance is > 100km. The 'cheap' measure is inaccurate over such
#> large distances, you'd likely be better using a different 'measure'.
#> [1] 10 10
```

All outputs are distances in metres, calculated with a variety of
spherical and elliptical distance measures. Distance measures currently
implemented are Haversine, Vincenty (spherical and elliptical)), the
very fast [mapbox cheap
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

### Benchmarks of geodesic accuracy

The `geodist_benchmark()` function - the only other function provided by
the `geodist` package - compares the accuracy of the different metrics
to the nanometre-accuracy standard of [Karney
(2013)](https://link.springer.com/content/pdf/10.1007/s00190-012-0578-z.pdf).

``` r
geodist_benchmark (lat = 30, d = 1000)
#>            haversine    vincenty       cheap
#> absolute 0.733504364 0.733504364 0.546134467
#> relative 0.002072017 0.002072017 0.001601954
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
rbenchmark::benchmark (
    replications = 10, order = "test",
    d1 <- geodist (x, measure = "cheap"),
    d2 <- geodist (x, measure = "haversine"),
    d3 <- geodist (x, measure = "vincenty"),
    d4 <- geodist (x, measure = "geodesic")
) [, 1:4]
#>                                      test replications elapsed relative
#> 1     d1 <- geodist(x, measure = "cheap")           10   0.068    1.000
#> 2 d2 <- geodist(x, measure = "haversine")           10   0.139    2.044
#> 3  d3 <- geodist(x, measure = "vincenty")           10   0.229    3.368
#> 4  d4 <- geodist(x, measure = "geodesic")           10   3.315   48.750
```

Geodesic distance calculation is available in the [`sf`
package](https://cran.r-project.org/package=sf). Comparing computation
speeds requires conversion of sets of numeric lon-lat points to `sf`
form with the following code:

``` r
require (magrittr)
x_to_sf <- function (x) {
    sapply (seq (nrow (x)), function (i) {
        sf::st_point (x [i, ]) %>%
            sf::st_sfc ()
    }) %>%
        sf::st_sfc (crs = 4326)
}
```

``` r
n <- 1e2
x <- cbind (-180 + 360 * runif (n), -90 + 180 * runif (n))
colnames (x) <- c ("x", "y")
xsf <- x_to_sf (x)
#> Warning in CPL_gdal_init(): GDAL Error 1: libpodofo.so.0.9.8: cannot open
#> shared object file: No such file or directory

#> Warning in CPL_gdal_init(): GDAL Error 1: libpodofo.so.0.9.8: cannot open
#> shared object file: No such file or directory

#> Warning in CPL_gdal_init(): GDAL Error 1: libpodofo.so.0.9.8: cannot open
#> shared object file: No such file or directory

#> Warning in CPL_gdal_init(): GDAL Error 1: libpodofo.so.0.9.8: cannot open
#> shared object file: No such file or directory
sf_dist <- function (xsf) sf::st_distance (xsf, xsf)
geo_dist <- function (x) geodist (x, measure = "geodesic")
rbenchmark::benchmark (
    replications = 10, order = "test",
    sf_dist (xsf),
    geo_dist (x)
) [, 1:4]
#>           test replications elapsed relative
#> 2  geo_dist(x)           10   0.070    1.000
#> 1 sf_dist(xsf)           10   0.144    2.057
```

Confirm that the two give almost identical results:

``` r
ds <- matrix (as.numeric (sf_dist (xsf)), nrow = length (xsf))
dg <- geodist (x, measure = "geodesic")
formatC (max (abs (ds - dg)), format = "e")
#> [1] "3.7676e+04"
```

All results are in metres, so the two differ by only around 10
nanometres.

The [`geosphere` package](https://cran.r-project.org/package=geosphere)
also offers sequential calculation which is benchmarked with the
following code:

``` r
fgeodist <- function () geodist (x, measure = "vincenty", sequential = TRUE)
fgeosph <- function () geosphere::distVincentySphere (x)
rbenchmark::benchmark (
    replications = 10, order = "test",
    fgeodist (),
    fgeosph ()
) [, 1:4]
#> The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
#> which was just loaded, will retire in October 2023.
#> Please refer to R-spatial evolution reports for details, especially
#> https://r-spatial.org/r/2023/05/15/evolution4.html.
#> It may be desirable to make the sf package available;
#> package maintainers should consider adding sf to Suggests:.
#> The sp package is now running under evolution status 2
#>      (status 2 uses the sf package in place of rgdal)
#>         test replications elapsed relative
#> 1 fgeodist()           10   0.017    1.000
#> 2  fgeosph()           10   0.032    1.882
```

`geodist` is thus around 3 times faster than `sf` for highly accurate
geodesic distance calculations, and around twice as fast as `geosphere`
for calculation of sequential distances.

### Test Results

``` r
require (devtools)
require (testthat)
```

``` r
date ()
#> [1] "Mon Aug 21 11:34:48 2023"
devtools::test ("tests/")
#> ℹ Testing geodist
#> Starting 2 test processes
#> ✔ | F W S  OK | Context
#> ⠋ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ] Starting up...                                                                                                                                                                                                                         ⠙ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ] Starting up...                                                                                                                                                                                                                         ⠹ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ] Starting up...                                                                                                                                                                                                                         ⠸ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ] Starting up...                                                                                                                                                                                                                         ⠼ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ] Starting up...                                                                                                                                                                                                                         ⠴ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ] Starting up...                                                                                                                                                                                                                         ⠦ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ] Starting up...                                                                                                                                                                                                                         ✔ |         4 | geodist-min                                                                                                                                                                                                                                                    
#> ⠧ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ] @ geodist                                                                                                                                                                                                                             ⠇ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ] @ geodist                                                                                                                                                                                                                             ✔ |        52 | geodist                                                                                                                                                                                                                                                        
#> ⠏ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 87 ] @ georange                                                                                                                                                                                                                            ⠋ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 87 ] @ georange                                                                                                                                                                                                                            ✔ |        37 | georange                                                                                                                                                                                                                                                       
#> ⠙ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 97 ] @ input-format                                                                                                                                                                                                                        ⠹ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 98 ] @ input-format                                                                                                                                                                                                                        ✔ |        18 | input-format                                                                                                                                                                                                                                                   
#> ⠸ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 111 ] @ measures                                                                                                                                                                                                                           ⠼ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 111 ] @ measures                                                                                                                                                                                                                           ✔ |        18 | measures                                                                                                                                                                                                                                                       
#> ⠴ [ FAIL 0 | WARN 0 | SKIP 0 | PASS 129 ] Starting up...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
#> ══ Results ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
#> Duration: 1.5 s
#> 
#> [ FAIL 0 | WARN 0 | SKIP 0 | PASS 129 ]
```

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the
[`allcontributors`
package](https://github.com/ropenscilabs/allcontributors) following the
[all-contributors](https://allcontributors.org) specification.
Contributions of any kind are welcome!

### Code

<table>
<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/commits?author=mpadge">mpadge</a>
</td>
<td align="center">
<a href="https://github.com/daniellemccool">
<img src="https://avatars.githubusercontent.com/u/5112209?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/commits?author=daniellemccool">daniellemccool</a>
</td>
</tr>
</table>

### Issues

<table>
<tr>
<td align="center">
<a href="https://github.com/mdsumner">
<img src="https://avatars.githubusercontent.com/u/4107631?u=c7a3627c592123651d51d002f421c2bd00be172f&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Amdsumner">mdsumner</a>
</td>
<td align="center">
<a href="https://github.com/edzer">
<img src="https://avatars.githubusercontent.com/u/520851?u=9bc892c3523be428dc211f2ccbcf04e8e0e564ff&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Aedzer">edzer</a>
</td>
<td align="center">
<a href="https://github.com/njtierney">
<img src="https://avatars.githubusercontent.com/u/6488485?u=3eacd57f61342d1c3cecd5c8ac741b1c4897e1de&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Anjtierney">njtierney</a>
</td>
<td align="center">
<a href="https://github.com/mkuehn10">
<img src="https://avatars.githubusercontent.com/u/4229651?u=b46361a6063e22f9d902672cb7e776918a1c6951&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Amkuehn10">mkuehn10</a>
</td>
<td align="center">
<a href="https://github.com/asardaes">
<img src="https://avatars.githubusercontent.com/u/7768461?u=fb573498b515f9bfcaeba8e256852060f6304d0b&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Aasardaes">asardaes</a>
</td>
<td align="center">
<a href="https://github.com/marcosci">
<img src="https://avatars.githubusercontent.com/u/10864574?u=e5b7e55e122646f47174a9e621ebc91fff177d9b&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Amarcosci">marcosci</a>
</td>
<td align="center">
<a href="https://github.com/mem48">
<img src="https://avatars.githubusercontent.com/u/15819577?u=0c128db4e7567656c23e83e4314111fcea424526&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Amem48">mem48</a>
</td>
</tr>
<tr>
<td align="center">
<a href="https://github.com/dcooley">
<img src="https://avatars.githubusercontent.com/u/8093396?u=2c8d9162f246d90d433034d212b29a19e0f245c1&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Adcooley">dcooley</a>
</td>
<td align="center">
<a href="https://github.com/Robinlovelace">
<img src="https://avatars.githubusercontent.com/u/1825120?u=461318c239e721dc40668e4b0ce6cc47731328ac&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3ARobinlovelace">Robinlovelace</a>
</td>
<td align="center">
<a href="https://github.com/espinielli">
<img src="https://avatars.githubusercontent.com/u/891692?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3Aespinielli">espinielli</a>
</td>
<td align="center">
<a href="https://github.com/Maschette">
<img src="https://avatars.githubusercontent.com/u/14663215?u=93694159d02e924e6413bd067d7746f1d16d64c1&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/hypertidy/geodist/issues?q=is%3Aissue+author%3AMaschette">Maschette</a>
</td>
</tr>
</table>
<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
