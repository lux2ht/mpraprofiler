
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mpraprofiler

<!-- badges: start -->

<!-- badges: end -->

The goal of mpraprofiler is to create a simple tool to perform allelic
differential expression analysis for **massively parallel reporter
assays** data

## Installation

You can install the the development version version of mpraprofiler from
gitlab with. It is still too early for this package to be public. You
can use `getPass` package to pass the authorization for gitLab:

``` r
# install.packages("devtools")
# install.packages("getPass")
devtools::install_git("https://tfwebdev.research.cchmc.org/gitlab/lux2ht/mpraprofiler.git", ,credentials = git2r::cred_user_pass("lux2ht", getPass::getPass()))
```

## Example

please refer to **vignettes** for details
