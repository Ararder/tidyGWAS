# Infer what genome build a GWAS summary statistics file is on.

Infer what genome build a GWAS summary statistics file is on.

## Usage

``` r
infer_build(tbl, dbsnp_path, n_snps = 10000)
```

## Arguments

- tbl:

  a `data.frame` or
  [`character()`](https://rdrr.io/r/base/character.html) vector

- dbsnp_path:

  filepath to the dbSNP155 directory

- n_snps:

  number of snps to check CHR and POS for

## Value

either "37" or "38"

## Examples

``` r
if (FALSE) { # \dontrun{
genome_build <- infer_build(gwas_sumstats)
} # }
```
