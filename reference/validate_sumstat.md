# Validate statistics columns in a GWAS summary statistics file

Validate statistics columns in a GWAS summary statistics file

## Usage

``` r
validate_sumstat(tbl, remove_cols = c(""), filter_func, convert_p)
```

## Arguments

- tbl:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

- remove_cols:

  Columns that should not be validated

- filter_func:

  handles reporting and writing removed files to disk

- convert_p:

  What value should be used for when P-value has been rounded to 0?

## Value

a tbl

## Examples

``` r
if (FALSE) { # \dontrun{
validate_sumstat(sumstat, remove_cols = "EffectAllele", convert_p = 0)
} # }
```
