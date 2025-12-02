# ![Logo](pkgdown/favicon/synthos_logo.png)

[![CRAN](https://www.r-pkg.org/badges/version/mbg?color=ffcc00)](https://cran.r-project.org/package=mbg)
[![Total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/mbg?color=blue)](https://cran.r-project.org/package=mbg)
[![Build
status](https://github.com/open-AIMS/synthos/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/open-AIMS/synthos/actions/workflows/pkgdown.yaml)

**`synthos` is an R package to generate synthetic data.**

The `synthos` package provides a simple interface to generate synthetic
data for ecological communities. ADD MORE TEXT

The `synthos` package combines features from the
[`sf`](https://r-spatial.github.io/sf/) and
[`stars`](https://r-spatial.github.io/stars/) packages for spatial data
processing; and, [`R-INLA`](https://www.r-inla.org/) and \[`gstat`\]
(<https://r-spatial.github.io/gstat/>) for geostatistical models.

------------------------------------------------------------------------

## Using the package

**You can install the latest version of the synthos package:**

`remotes::install_github("open-AIMS/synthos@julie")`

Some core package functions rely on R-INLA, which is not available on
CRAN. If you do not already have the `INLA` package installed, you can
download it at (<https://www.r-inla.org/download-install>).

After installing and package and loading it using
[`library(synthos)`](https://open-aims.github.io/synthos/), you can
access the package vignette by running `help(mbg)`, or get documentation
for a specific function by running e.g. `help(MbgModelRunner)`.

------------------------------------------------------------------------
