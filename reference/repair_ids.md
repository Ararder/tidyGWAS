# Augment a data.frame with information from dbSNP

Augment a data.frame with information from dbSNP

## Usage

``` r
repair_ids(
  tbl,
  repair = c("rsid", "pos"),
  build = c("NA", "37", "38"),
  dbsnp_path
)
```

## Arguments

- tbl:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with columns CHR, POS, EffectAllele, OtherAllele or RSID,
  EffectAllele, OtherAllele2.

- repair:

  "rsid" to repair RSID, "pos" to repair CHR and POS

- build:

  used if repair = "rsid" to specify genome build

- dbsnp_path:

  path to the directory containing the dbSNP dataset

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
with added columns CHR, POS_37, POS_38 RSID, REF_37, ALT_37, REF_38,
ALT_38, no_dbsnp_entry, incompat_alleles

## Examples

``` r
if (FALSE) { # \dontrun{
repair_ids(gwas_sumstats, repair = "rsid", build = "38", dbsnp_path = "/path/to/dbsnp")
} # }
```
