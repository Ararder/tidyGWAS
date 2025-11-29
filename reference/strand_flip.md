# Strand flip alleles

Strand flip alleles

## Usage

``` r
strand_flip(tbl)
```

## Arguments

- tbl:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with columns `EffectAllele` and `OtherAllele`

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
with columns `EffectAllele` and `OtherAllele` flipped

## Examples

``` r
if (FALSE) { # \dontrun{
tbl <- strand_flip(tbl)
} # }
```
