# Estimate ancestry composition from allele frequencies

This function implements the `bigsnpr::snp_ancestry_summary()` written
by Florian Prive.

## Usage

``` r
ancestry_comp(tbl, dbsnp_path, min_cor = 0.4, sum_to_one = TRUE)
```

## Arguments

- tbl:

  a data.frame with atleast columns `RSID`, `EffectAllele`,
  `OtherAllele`, and `EAF`.

- dbsnp_path:

  filepath to the dbSNP reference data. Same path as for
  [`tidyGWAS()`](https://ararder.github.io/tidyGWAS/reference/tidyGWAS.md)

- min_cor:

  minimum correlation between predicted and observed allele frequencies.

- sum_to_one:

  Force ancestry coefficients to sum to 1?

## Value

a tibble with columns `est` and `ancestry`

## Details

The function and reference data is based on
<https://doi.org/10.1093/bioinformatics/btac348>

`ancestry_comp()` requires additional reference data. This data is only
available if you have downloaded the updated reference tidyGWAS data The
reference file `ancestry_data.parquet` is expected to exist in the same
dbSNP folder that is used with
[`tidyGWAS()`](https://ararder.github.io/tidyGWAS/reference/tidyGWAS.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ancestry_comp(tbl, "path_to_data")

} # }
```
