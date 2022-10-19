# v0.0.8

### Minor changes:

- Update geodesic source code
- Add codeberg.org/hypertidy/geodist remote
- Fix all clang warnings, including new header files with prototypes for all fns

# v0.0.7

### Minor changes:

- Improve matching of column names to lon/lat
- Fix error when submitting tibbles (#32)


# v0.0.6

### Minor changes:

- Mirrored source repository on gitlab, bitbucket, sourcehut, and github
- All functions now issue a message to use a different measure when default
  `measure = "cheap"` and maximum distance > 100km (see #26)
- C code made more robust (all input variables protected).

# v0.0.4

### Major changes:

- Add 'geodist_vec()' function to accept individual vectors

### Minor changes:

- bug fix for cheap distances between single points; thanks to @mem48

# v0.0.3

### Minor changes:

- bug fix in sequential calculation with cheap distances; thanks to @mem48

# v0.0.2

### Major changes:

- Add georange function to calculate extreme ranges only; thanks to @marcosci
- Add 'paired' option to 'geodist()' function to calculate pair-wise distances
  between 'x' and 'y' matrices

### Minor changes:

- Improve handling of column names


# v0.0.1

Initial CRAN release.
