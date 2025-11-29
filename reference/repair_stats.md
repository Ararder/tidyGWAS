# Repair statistics column in a GWAS summary statistics tibble

`repair_stats()` is a collection of functions that can be used to infer
missing columns in GWAS summary statistics. The functions are based on
functionality found online.

## Usage

``` r
repair_stats(
  tbl,
  dbsnp_path,
  impute_freq = c("None", "EUR", "AMR", "AFR", "SAS", "EAS"),
  impute_freq_file = NULL,
  impute_n = FALSE
)
```

## Arguments

- tbl:

  a `data.frame` or
  [`character()`](https://rdrr.io/r/base/character.html) vector

- dbsnp_path:

  filepath to the dbSNP155 directory

- impute_freq:

  one of c("None", "EUR", "AMR", "AFR", "SAS", "EAS"). If None, no
  imputation is done. Otherwise precomputed alleles frequence from
  1000KG, selected ancestry is used

- impute_freq_file:

  filepath to a .parquet file with custom allele frequencies. The file
  needs to be a tabular dataframe with columns RSID, EffectAllele,
  OtherAllele, EAF. EAF should correspond to the frequency of the
  EffectAllele.

- impute_n:

  Should N be imputed if it's missing?

## Value

a tibble

## Examples

``` r
if (FALSE) { # \dontrun{
updated <- repair_stats(my_gwas)
} # }
```
