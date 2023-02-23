
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rmexact

<!-- badges: start -->

[![NSF-2132247](https://img.shields.io/badge/NSF-2132247-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=2132247)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/gerardlab/rmexact/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gerardlab/rmexact/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/gerardlab/rmexact/branch/main/graph/badge.svg)](https://app.codecov.io/gh/gerardlab/rmexact?branch=main)
[![DOI](https://zenodo.org/badge/525824572.svg)](https://zenodo.org/badge/latestdoi/525824572)
<!-- badges: end -->

Exact tests for random mating in autopolyploids. The main functions are

- `tetexact()`: Exact test using sufficient statistics for random mating
  in autotetraploids when the genotypes are known.
- `rmslrt()`: Exact test using the split likelihood ratio approach of
  Wasserman et al. (2020) for random mating in autopolyploids when
  genotypes are known.
- `rmchisq()`: Chi-squared test for random mating in autopolyploids when
  genotypes are known.

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

## References

Wasserman, L., Ramdas, A., & Balakrishnan, S. (2020). Universal
inference. *Proceedings of the National Academy of Sciences*, 117(29),
16880-16890.
[doi:10.1073/pnas.1922664117](https://doi.org/10.1073/pnas.1922664117)
