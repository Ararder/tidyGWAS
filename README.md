
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyGWAS

<!-- badges: start -->
<!-- badges: end -->

The goal of tidyGWAS is to provide an easy-to use R package written in
the typical tidyverse/dplyr syntax to make cleaning and collection of
sumstats simple, and to stop analysts for falling to common error. In
addition we provide helpers to standardise GWAS data: Reparation of
stats columns where possible (missing B/SE/P), as well as by default
providing coordinates on both GRCh37 and GRCh38.

An early design goal of tidyGWAS was to keep as many rows as possible.
Our default options are therefore to keep multi-allelic SNPs or
ambigious SNPs, which we suggest to handle depending on what downstream
analysis tool you wish to use.

## Installation

You can install the development version of tidyGWAS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Ararder/tidyGWAS")
```

## Example

see vignette for “tidyGWAS”
