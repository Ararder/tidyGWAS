# Query a specific region of interest for a using a gwas catalog study_id

Query a specific region of interest for a using a gwas catalog study_id

## Usage

``` r
from_gwas_catalog_region(study_id, chr, start, end)
```

## Arguments

- study_id:

  a study accession ID, e.g. "GCST000001"

- chr:

  chromosome number, e.g. "1" - not "chr1"

- start:

  base pair start position, e.g. 1000000

- end:

  base pair end position, e.g. 2000000

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

## Examples

``` r
if (FALSE) { # \dontrun{
#' get_gwas_catalog_region("GCST000001", "1", 1000000, 2000000)
} # }
```
