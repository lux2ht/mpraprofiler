
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mpraprofiler

<!-- badges: start -->
<!-- badges: end -->

The goal of mpraprofiler is to create a simple tool to perform allelic
differential expression analysis for massively parallel reporter assay
(MPRA) data.

## Installation

You can install mpraprofiler from
[GitHub](https://github.com/WeirauchLab/mpraprofiler) with:

``` r
install.packages("devtools")  # if not already installed
devtools::install_github("WeirauchLab/mpraprofiler")
```

To build the vignette(s), install some additional packages

``` r
install.packages("tidyverse")
install.packages("readxl")

# You'll need Bioconductor to get DESeq2
install.packages("BiocManager")
BiocManager::install('DESeq2')
```

…and add `build_vignettes = TRUE` to the install command:

``` r
devtools::install_github("WeirauchLab/mpraprofiler", build_vignettes = TRUE,
                         force = TRUE)
```

## Example Usage

Please refer to
[`vignettes/sample_analysis.md`](vignettes/sample_analysis.md) for
details.

## License

GNU General Public License v3. See [`LICENSE.md`](LICENSE.md)
