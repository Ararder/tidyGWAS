# Meta analyse a chromosome, one at a time!

Meta analyse a chromosome, one at a time!

## Usage

``` r
by_chrom(ds, chrom, ref)
```

## Arguments

- ds:

  an
  [`arrow::open_dataset()`](https://arrow.apache.org/docs/r/reference/open_dataset.html)
  object

- chrom:

  chromosome

- ref:

  reference allele

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

## Examples

``` r
if (FALSE) { # \dontrun{
by_chrom(ds, chrom = 1, ref = "REF_38")
} # }
```
