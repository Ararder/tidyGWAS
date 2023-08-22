
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
automatic standardisation and quality control of summary statistics. See
the package [webpage](http://arvidharder.com/tidyGWAS/) to get started

1.  Detection of duplicated rows (based on RSID_REF_ALT or
    CHR_POS_REF_ALT)

2.  Standardized column names

3.  Automatic updating of “merged” RSIDs

4.  Detection and optional removal of deletions/insertions (“indels”)

5.  Automatic detection and conversion of CHR:POS or CHR:POS:A1:A2 to
    CHR, POS, RSID, A1 and A2.

6.  Standardization of CHR values (“23” -\> “X”, “chr1” -\> “1”)

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
