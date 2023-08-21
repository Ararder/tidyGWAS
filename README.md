
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyGWAS

<!-- badges: start -->
<!-- badges: end -->

# Standardized and automatic cleaning and harmonization of GWAS summary statistics

The typical process of starting a project that relies on GWAS summary
statistics looks like this:

1.  Read through data access statements to find the downloadable link
2.  Dig through readme files to identify genome build and what exactly
    column names mean
3.  Manually loading the sumstats into R/Python/Bash, taking a look at
    it to identify any obvious errors
4.  For each sumstat, format it into whatever format is requested by the
    downstream tool you are using
5.  Rinse and repeat for the next project

To solve this issue we developed a R package “tidyGWAS”, that does
automatic standardisation and quality control of summary statistics.

1.  Detection of duplicated rows (based on RSID_REF_ALT or
    CHR_POS_REF_ALT)

2.  Standardized column names

3.  Automatic updating of “merged” RSIDs

4.  Detection and optional removal of deletions/insertions (“indels”)

5.  Automatic detection and conversion of CHR:POS or CHR:POS:A1:A2 to
    CHR, POS, RSID, A1 and A2.

6.  Standardization of CHR values (23 -\> “X”, chr1 -\> “1”)

7.  Validation of standard GWAS columns, B, SE, P, N, FREQ, Z, CaseN,
    ControlN,A1, A2

    1.  Extremely small pvalues are by default converted to
        2.225074e-308 (minimum pvalue in R)

8.  Imputation of missing columns: RSID from CHR:POS or CHR:POS from
    RSID. Any of B,SE, P, Z if missing and possible

9.  Validation of CHR:POS:RSID by matching with dbSNP v.155

10. Cleaned sumstats are provided with coordinates on both GRCh37 and
    GRCh38, with TRUE/FALSE flags for indels and variants that are
    multi-allelic in the dataset

# Getting started

tidyGWAS dependencies should be easily installed on most systems, with
the arrow dependecy being the most likely to cause trouble. In the case
that installation of arrow is causing trouble, see
[here](https://arrow.apache.org/docs/r/articles/install.html).

``` r
# install.packages("devtools")
devtools::install_github("ararder/tidyGWAS")
# or
# install.packages("remotes")
remotes::install_github("ararder/tidyGWAS")
```

tidyGWAS uses a version dbSNP 155 converted to the apache arrow parquet
files. You can download the dbSNP reference file from inside R, or by
manually navigating to this
[file](https://drive.google.com/file/d/1LmgfmTQaWwJaFpBHcRQIY_kwe5iN7Pj6/view?usp=share_link)

``` r
# You can download the file from inside R like this:
library(googledrive)
id <- googledrive::as_id("1LmgfmTQaWwJaFpBHcRQIY_kwe5iN7Pj6")

##### EDIT THIS:
filepath_to_store_dir <- ""
##### ---------------------

drive_download(id, path = filepath_to_store_dir)
```

The file needs to be untarred to use

``` bash
cd $filepath_to_store_dir
tar -xf dbSNP155
```

## Example use

``` r
library(tidyGWAS)
dbsnp_file <-  "filepath/to/untarred/dnsnp155"
# input a filepath, or a (data.frame / tibble / data.table), of summary statistics

# by default, tidyGWAS returns the cleaned sumstats inside R
# AND writes out the file in a temp directory
cleaned <- tidyGWAS(
  tbl = tidyGWAS::example_file,
  dbsnp_files = "filepath/to/dbSNP155"
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
