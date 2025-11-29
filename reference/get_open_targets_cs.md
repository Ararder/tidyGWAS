# Query Open targets for all credible sets containing the variant

Query Open targets for all credible sets containing the variant

## Usage

``` r
get_open_targets_cs(
  variant_id,
  page_size = 500,
  api_url = "https://api.platform.opentargets.org/api/v4/graphql"
)
```

## Arguments

- variant_id:

  in the open targets format chr_pos_ref_alt: 5_100_101_A_C

- page_size:

  number of results per page, default 500

- api_url:

  URL of the Open Targets GraphQL API, default

## Value

a [dplyr::tibble](https://dplyr.tidyverse.org/reference/reexports.html)
with the following columns:

- `CHR`: chromosome

- `POS`: position

- `P`: p-value

- `REF`: reference allele

- `ALT`: alternate allele

- `description`: trait description

- `study_id`: study ID

## Examples

``` r
if (FALSE) { # \dontrun{
get_open_targets_cs("7_140459051_C_G")
} # }
```
