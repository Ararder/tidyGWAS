# Detect entries that are not valid rsID's in GWAS summary statistics

Detect entries that are not valid rsID's in GWAS summary statistics

## Usage

``` r
flag_invalid_rsid(tbl, regex = "^[rR][sS]?\\d{1,10}$")
```

## Arguments

- tbl:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with column `RSID`.

- regex:

  regex used to detect non-RSIDs

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
with column `invalid_rsid`

## Examples

``` r
if (FALSE) { # \dontrun{
flag_invalid_rsid(tbl) |>
dplyr::filter(invalid_rsid)
} # }
```
