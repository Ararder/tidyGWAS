# Perform meta-analysis of GWAS summary statistics datasets cleaned by tidyGWAS This function is depreciated. Use [`meta_analyse()`](https://ararder.github.io/tidyGWAS/reference/meta_analyse.md)

Perform meta-analysis of GWAS summary statistics datasets cleaned by
tidyGWAS This function is depreciated. Use
[`meta_analyse()`](https://ararder.github.io/tidyGWAS/reference/meta_analyse.md)

## Usage

``` r
deprec_meta_analyze(
  dset,
  by = c("CHR", "POS_37", "RSID", "EffectAllele", "OtherAllele"),
  ref = c("REF_37", "REF_38")
)
```

## Arguments

- dset:

  an
  [`arrow::open_dataset()`](https://arrow.apache.org/docs/r/reference/open_dataset.html)
  object

- by:

  a character vector of column names to group by. Default is c("CHR",
  "POS", "RSID", "EffectAllele", "OtherAllele")

- ref:

  either "REF_37" or "REF_38", depending on which column you want to use
  to standardize reference allele

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

## Examples

``` r
if (FALSE) { # \dontrun{
dset <- arrow::open_dataset("path_to/sumstats/")
res <- meta_analyze(dset)
} # }
```
