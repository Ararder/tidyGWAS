# Validate format of the RSID column in a GWAS summary statistics file

Validate format of the RSID column in a GWAS summary statistics file

## Usage

``` r
validate_rsid(tbl, filepath)
```

## Arguments

- tbl:

  a `data.frame` or
  [`character()`](https://rdrr.io/r/base/character.html) vector

- filepath:

  filepath to write out removed rows

## Value

a tbl

## Examples

``` r
if (FALSE) { # \dontrun{
validate_rsid(sumstat, "~/invalid_rsid.parquet")
} # }
```
