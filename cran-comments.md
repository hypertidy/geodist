# CRAN notes for geodist_0.0.8 submission

The current version (0.0.7) generates Warnings on all CRAN Linux machines: "a function declaration without a prototype is deprecated in all versions of C". The function declaration referred to is, however, *static*, so does not require a prototype. Moreover, I can not replicate that warning with Clang v14.0.6 using "-Wstrict-prototype".

Other than that, this submission generates no errors, warning or notes on the following R versions and operating systems:

* Ubuntu: R-oldrelease, R-release, R-devel
- OSX: R-release
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
