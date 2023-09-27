
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyGWAS

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/Ararder/tidyGWAS/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Ararder/tidyGWAS?branch=main)
<!-- badges: end -->

Genome-wide summary statistics are becoming a staple in many different
genetics and genomics analysis pipelines. Often, the specific filters
suggested for pipelines can be different, requiring each pipeline to
have a step where summary statistics are “munged”.

tidyGWAS aims to provide a standardized format *before* before any
pipeline specific munging is done. With that in mind, tidyGWAS is
conservative in removing rows, and by default keeps both indels and
multi-allelic variants.

`tidyGWAS` does the following:

1.  Detection of duplicated rows (based on RSID_REF_ALT or
    CHR_POS_REF_ALT)

2.  Standardized column names

3.  Automatic updating of
    [merged](https://www.ncbi.nlm.nih.gov/books/NBK573473/) RSIDs

4.  Detection and optional removal of deletions/insertions (“indels”)

5.  Detection of non rsID values in RSID column, and automatic parsing
    of the common CHR:POS or CHR:POS:REF:ALT format

6.  Standardization of CHR values (ex: “23” -\> “X”, “chr1” -\> “1”)

7.  Validation of standard GWAS columns, B, SE, P, N, FREQ, Z, CaseN,
    ControlN, A1, A2

    1.  Extremely small pvalues are by default converted to
        2.225074e-308 (minimum pvalue in R)

8.  Imputation of missing columns: RSID from CHR:POS or CHR:POS from
    RSID. Any of B,SE, P, Z if missing and possible

9.  Validation of CHR:POS:RSID by matching with dbSNP v.155

10. Cleaned sumstats are provided with coordinates on both GRCh37 and
    GRCh38, with TRUE/FALSE flags for indels and variants that are
    multi-allelic in the dataset

From working with standardized GWAS formats, we’ve found that having
both GRCh37 and GRCh38 coordinates, and standardized column names
significantly speeds up downstream analysis.

The computationally intensive part of aligning summary statistics with
dbSNP 155 (\> 940 million rows) for both GRCh37 and GRCh38 (in total 1.8
billion rows) is implemented using the [Apache Arrow
R](https://arrow.apache.org/docs/r/) implementation, allowing for the
full function to run in \<5 minutes, using less than 16gb, with ~7
million rows on a Macbook Pro M2.
