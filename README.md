
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BREADR

<!-- badges: start -->

[![R-CMD-check](https://github.com/jonotuke/BREADR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jonotuke/BREADR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of BREADR is to provide an easy-to-use method for estimating
degrees of relatedness (up to the second degree) for extremely
low-coverage data. BREADR also allows users to quantify and visualise
the level of confidence in the estimated degrees of relatedness.

The method requires Eigenstrat files (an ind, geno and snp file) to
begin, allowing the user to account for DNA deamination when genotyping
their data.

Further information can be found at the [BREADR
website](https://jonotuke.github.io/BREADR/).

## Installation

To install, you can use the usual

``` r
install.packages("BREADR")
```

You can install the development version of BREADR from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("jonotuke/BREADR")
```

Also to ensure that BREADR install correctly, we suggest installing the
following packages

- Matrix,
- data.table,
- dplyr,
- forcats,
- ggplot2,
- ggpubr,
- grDevices,
- magrittr,
- MASS,
- matrixStats,
- purrr,
- readr,
- stringr and
- tibble.

## Contributing to BREADR

Please note that `BREADR` is work in progress! If you are interested in
the project and want to know more please feel free to contact Jono Tuke
(<simon.tuke@adelaide.edu.au>). If you find a bug or would like to see
new or improved features, please open an issue on the
[GitHub](https://github.com/jonotuke/BREADR) repository.
