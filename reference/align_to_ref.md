# Align EffectAllele to always be the reference genome allele

`align_to_ref()` will:

1.  Filter any variants where EffectAllele or OtherAllele is not the
    reference genome allele

2.  Flip EffectAllele such that it is always the reference genome allele

3.  If OtherAllele is the reference genome allele, direction of B,Z and
    EAF are flipped.

## Usage

``` r
align_to_ref(dset, ref = c("REF_38", "REF_37"))
```

## Arguments

- dset:

  object created by
  [`arrow::open_dataset()`](https://arrow.apache.org/docs/r/reference/open_dataset.html)
  or
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

- ref:

  Reference genome allele to align to. Varies in ~0.5% of locations

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
or arrow query depending on whether dset is a tibble or arrow dataset

## Examples

``` r
if (FALSE) { # \dontrun{
align_to_ref(tidygwas_df, ref = "REF_38")
} # }
```
