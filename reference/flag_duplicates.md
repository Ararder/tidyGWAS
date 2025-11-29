# Find all rows which are part of a set of duplicated rows

Many duplication tools such as
[`base::duplicated()`](https://rdrr.io/r/base/duplicated.html) or
[`dplyr::distinct()`](https://dplyr.tidyverse.org/reference/distinct.html)
identify rows which are duplications. It is often useful to see ALL rows
which are part of the duplication set, and not just the second row.

creates new column: `dup_rsid` or `dup_chr_pos`, a T/F flag.
Specifically, flags both rows in a duplication pair, and not just first
or last duplicate row, making it easy to work with all rows that are
part of a duplication

## Usage

``` r
flag_duplicates(
  tbl,
  column = c("rsid", "chr_pos", "chr_pos_ref_alt", "rsid_ref_alt")
)
```

## Arguments

- tbl:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

- column:

  Which columns should be used to form a unique ID?

## Value

a tibble with a new column marking duplicates

## Examples

``` r
if (FALSE) { # \dontrun{

# will tag multi-allelics as duplications
flag_duplicates(tbl, column = "rsid")
flag_duplicates(tbl, column = "chr_pos")
# if you are interested in rows that are variant duplications
flag_duplicates(tbl, column = "rsid_ref_alt")
flag_duplicates(tbl, column = "chr_pos_ref_alt")

} # }
```
