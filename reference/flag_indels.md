# Detect "indels" in GWAS summary statistics

Detect "indels" in GWAS summary statistics

## Usage

``` r
flag_indels(tbl)
```

## Arguments

- tbl:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with columns `EffectAllele` and `OtherAllele`

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
with a TRUE/FALSE column `indel` added, where indel == TRUE corresponds
to a row marked as an indel.

## Examples

``` r
if (FALSE) { # \dontrun{
all_indels <-
  flag_indels(tbl) |>
  dplyr::filter(indels)
} # }
```
