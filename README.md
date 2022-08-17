
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rmexact

<!-- badges: start -->

[![NSF-2132247](https://img.shields.io/badge/NSF-2132247-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=2132247)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/dcgerard/rmexact/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dcgerard/rmexact/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Exact tests for random mating in autopolyploids. We allow for genotype
uncertainty through the use of genotype likelihoods and fuzzy
*p*-values. The main functions are

-   `tetexact()`: Exact test for random mating in autotetraploids when
    the genotypes are known.
-   `tetgp()`: Exact test for random mating in autotetraploids using
    genotype posteriors.
-   `vpv()`: Calculate the vacuumed *p*-value.

## Installation

You can install the development version of `{rmexact}` using the
`{remotes}` package:

``` r
# install.packages("remotes")
remotes::install_github("dcgerard/rmexact")
```

## Code of Conduct

Please note that the rmexact project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
