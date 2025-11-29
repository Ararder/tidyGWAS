# Execute validation and quality control of GWAS summmary statistics

`tidyGWAS()` performs a set of validations on input colummns, repairs
missing columns, and can add missing CHR/POS or RSID. In addition, CHR
and POS is standardised to GRCh38, with coordinates on GRCh37 added in
as well.

Briefly, `tidyGWAS()` updates RSID if possible using the
[refsnp-merged](https://ftp.ncbi.nih.gov/snp/latest_release/JSON/) file
from dbSNP. Each inputed column is then validated and coerced to the
correct type.

If statistis such as `P`, `B` are missing, `tidyGWAS()` will attempt to
impute them if possible using
[`repair_stats()`](https://ararder.github.io/tidyGWAS/reference/repair_stats.md)

Standard column names are assumed, BEFORE inputting into the function.
This is a deliberate decision as automatic parsing of some important
column names can be ambiguous For example, in some sumstats, A1 referes
to effect allele, while other formats use A1 as non-effect allele.

## Usage

``` r
tidyGWAS(
  tbl,
  dbsnp_path,
  ...,
  column_names = NULL,
  output_format = c("hivestyle", "parquet", "csv"),
  output_dir = tempfile(),
  CaseN = NULL,
  ControlN = NULL,
  N = NULL,
  est_ancestry = FALSE,
  impute_freq = c("None", "EUR", "AMR", "AFR", "SAS", "EAS"),
  impute_freq_file = NULL,
  impute_n = FALSE,
  min_EAF = NULL,
  flag_discrep_freq = c("None", "EUR", "AMR", "AFR", "SAS", "EAS"),
  allow_duplications = FALSE,
  build = c("NA", "37", "38"),
  default_build = c("37", "38"),
  indel_strategy = c("keep", "qc", "remove"),
  convert_p = 2.225074e-308,
  repair_cols = TRUE,
  logfile = FALSE
)
```

## Arguments

- tbl:

  a `data.frame` or
  [`character()`](https://rdrr.io/r/base/character.html) vector

- dbsnp_path:

  filepath to the dbSNP155 directory

- ...:

  pass additional arguments to
  [`arrow::read_delim_arrow()`](https://arrow.apache.org/docs/r/reference/read_delim_arrow.html),
  if tbl is a filepath.

- column_names:

  a named list of column names: `list(RSID = "SNP", POS = "BP")`

- output_format:

  How should the finished cleaned file be saved?

  - 'csv' corresponds to
    [`arrow::write_csv_arrow()`](https://arrow.apache.org/docs/r/reference/write_csv_arrow.html)

  - 'parquet' corresponds to
    [`arrow::write_parquet()`](https://arrow.apache.org/docs/r/reference/write_parquet.html)

  - 'hivestyle' corresponds to
    [`arrow::write_dataset()`](https://arrow.apache.org/docs/r/reference/write_dataset.html)
    split by `CHR`

- output_dir:

  filepath to a folder where tidyGWAS output will be stored. The folder
  should not yet exist. Note that the default argument is
  [`tempfile()`](https://rdrr.io/r/base/tempfile.html), meaning that
  tidyGWAS output will not be saved by default over R sessions.

- CaseN:

  manually input number of cases

- ControlN:

  manually input number of controls

- N:

  manually input sample size

- est_ancestry:

  Should
  [`ancestry_comp()`](https://ararder.github.io/tidyGWAS/reference/ancestry_comp.md)
  be run at the end of tidyGWAS to estimate the ancestry of the summary
  statistics?

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

- min_EAF:

  Apply a filter on allele frequency prior to applying the algorithm.
  Useful to speed up cleaning of very large files

- flag_discrep_freq:

  Should variants with allele frequency discrepancies be flagged?

- allow_duplications:

  Should duplicated variants be allowed? Useful if the munged sumstats
  are QTL sumstats

- build:

  If you are sure of what genome build ('37' or '38'), can be used to
  skip
  [`infer_build()`](https://ararder.github.io/tidyGWAS/reference/infer_build.md)
  and speed up computation

- default_build:

  If only RSID exists, the build cannot be inferred. Nonetheless,
  tidyGWAS applies a filter on incompatible alleles with GRCh37/38. In
  such a case, tidyGWAS needs to decide on which reference genome to
  compare alleles with.

- indel_strategy:

  Should indels be kept or removed?

- convert_p:

  What value should be used for when P-value has been rounded to 0?

- repair_cols:

  Should any missing statistical columns be repaired if possible? calls
  [`repair_stats()`](https://ararder.github.io/tidyGWAS/reference/repair_stats.md)
  if TRUE

- logfile:

  Should messages be redirected to a logfile?

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

## Examples

``` r
if (FALSE) { # \dontrun{
tidyGWAS(
  tbl = "path/to/GWAS_trait_X_.tsv.gz", logfile = TRUE,
  output_dir = "/store/GWAS/tidyGWAS/trait_X"
  )
} # }
```
