
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyGWAS

<!-- badges: start -->
<!-- badges: end -->

The goal of tidyGWAS is to provide an easy-to use R package written in
the typical tidyverse/dplyr syntax to make cleaning and collection of
sumstats simple, and to stop analysts for falling to common error. In
addition we provide helpers to standardise GWAS data: Reparation of
stats columns where possible (missing B/SE/P), as well as by default
providing coordinates on both GRCh37 and GRCh38.

An early design goal of tidyGWAS was to keep as many rows as possible.
Our default options are therefore to keep multi-allelic SNPs or
ambigious SNPs, which we suggest to handle depending on what downstream
analysis tool you wish to use.

## Installation

You can install the development version of tidyGWAS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Ararder/tidyGWAS")
```

``` bash
# and then download the following file
# from google drive: https://drive.google.com/file/d/1hg6jxQUj6UmEdLIH46md7eJcgfVlzmNn/view?usp=share_link
# which contains dbSNP 155 for GRCh37 and GRCH38 in parquet format
```

## Example

see vignette

``` r
library(tidyGWAS)
dbsnp_file <-  "filepath/to/untarred/dnsnp155"
raw_sumstats <- "filepath/to/raw_sumstats.csv" 
cleaned <- tidyGWAS(
  tbl = raw_sumstats,
  dbsnp_files = "filepath/to/raw_sumstats.csv",
  output_format = "parquet",
  name = "my_GWAS"
)
```

## What does tidyGWAS do?

2.  Remove rows with NAs
3.  Handles duplicates across CHR:POS:REF:ALT or RSID:REF:ALT. Note that
    this will NOT remove multi-allelics
4.  Can detect and remove indels. Control using `keep_indels`
5.  Update only RSIDs that have been merged into newer ones, using
    RsMergeArch
6.  Check that RSID follows correct format.
7.  Check for CHR:POS or CHR:POS:REF:ALT in invalid RSIDs, and automatic
    parsing of this format
8.  Reparation of of RSID for rows with invalid RSID, if format is
    CHR:POS or CHR:POS:REF:ALT
9.  Validation CHR and POS, automatic conversion of common CHR formats
    such as ‘X’, ‘M’, ‘MT’ or ‘23’, “chr6” into “6”.
10. Validation of statistics columns: B, Z, SE, P, N, EAF. See
    `stats_validation`
11. Validation of CHR, POS and RSID columns. Identify rows with
    inconsistent CHR, POS RSID combination.
12. Repair missing CHR/POS from RSID, or RSID from CHR:POS.
13. Remove rows where alleles dont match possible alleles in dbSNP
14. Remove RSIDs which are not in dbSNP
15. Add CHR and POS from both GRCh37 and GRCh38
