
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/xcatala/methclust.svg?branch=master)](https://travis-ci.com/xcatala/methclust)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/xcatala/methclust?branch=master&svg=true)](https://ci.appveyor.com/project/xcatala/methclust)
[![R build
status](https://github.com/xcatala/methclust/workflows/R-CMD-check/badge.svg)](https://github.com/xcatala/methclust/actions)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/methclust.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/methclust)
[![Codecov test
coverage](https://codecov.io/gh/xcatala/methclust/branch/master/graph/badge.svg)](https://codecov.io/gh/xcatala/methclust?branch=master)
<!-- badges: end -->

# methclust

The goal of methclust is to …

## Installation

You can install the released version of methclust from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("methclust")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xcatala/methclust")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(methclust)
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
