# Improved meta-analysis using tidyGWAS:ed files

`meta_analyse()` will:

- Flip the effect allele to be the reference allele (using either GRCh37
  or GRCh38), control with `ref`

- variant id is constructed using RSID:EffectAllele:OtherAllele (Which
  is now RSID:REF:ALT)

- EAF (allele frequency) and INFO (imputation quality) are weighted by
  sample size, if present

- CaseN, N, ControlN and EffectiveN are all summed and carried forward

## Usage

``` r
meta_analyse(
  ds,
  chromosomes = c(1:22),
  min_EAF = NULL,
  ref = c("REF_38", "REF_37"),
  safe_mode = FALSE
)
```

## Arguments

- ds:

  a dataset object, see
  [`arrow::open_dataset()`](https://arrow.apache.org/docs/r/reference/open_dataset.html)

- chromosomes:

  Which chrosomes to apply meta-analysis across. Default is autosomes:
  1:22

- min_EAF:

  Filter on minimal EAF (effect allele frequency) value. (0 - 0.5)

- ref:

  Reference genome version to use for reference allele

- safe_mode:

  Apply additional filters to ensure no missing/inf/NaN across key
  columns? This is set to FALSE by default because these filters should
  already be guaranteed with the use of
  [`tidyGWAS()`](https://ararder.github.io/tidyGWAS/reference/tidyGWAS.md)

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

## Examples

``` r
if (FALSE) { # \dontrun{
meta_analyse2(ds, min_EAF = 0.01)
} # }
```
