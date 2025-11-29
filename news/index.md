# Changelog

## tidyGWAS 1.0.0

This version comes with major changes to tidyGWAS.

1.  File-parsing is now done with data.table. No need to use the `delim`
    argument
2.  tidyGWAS now accepts four possible inputs to `tbl`

&nbsp;

1.  local filepath,
2.  in-memory data.frame
3.  a GWAS catalogs study id (will download and process)
4.  a URL to a downloadable file (will download and process)

&nbsp;

3.  Column name guessing is now done automatically, and in most cases
    you should not need to provide the column_names argument
4.  INDELS are now supported, but requires you to update the reference
    data, see the zenodo repository. Set indel_strategy =“qc” to apply
    the same QC strategy to indels
5.  tidyGWAS can now estimate ancestry proportions of the cleaned
    summary statistics by the function decribed in
    <https://academic.oup.com/bioinformatics/article/38/13/3477/6590645>

## tidyGWAS 0.9.9.2

Updates: 1. Added experimental support for automatic column name
guessing 2. Added experimental support for downloading and cleaning
files directly from GWAS catalog, using only study_id

## tidyGWAS 0.9.9

Added reference allele frequencies to the dbSNP155 reference data from
zenodo. Added the options to impute allele frequenecy from this
reference data.

## tidyGWAS 0.9.7

added functionality for QTL summary statistics by providing the option
to not remove duplicated variants.

in addition, added the option to provide a custom reference file for
allele frequency, that will impute the allele frequency from the
reference file if allele frequency is missing.

Added the option to impute N using SE and allele frequency.

## tidyGWAS 0.9.6

Added values that tidyGWAS detects as correct columns: 24 -\> “Y” 25 -\>
“XY” 26 -\> “MT”

## tidyGWAS 0.9.0

Major revision from earlier versions. The reference data format has been
updated:

1.  Compression switched from gzip to snappy (better speed, and less
    installation issues)
2.  GRCh37 and GRCh38 have been merged into a single file, requiring
    only one merge to get CHR, POS and RSID for both builds.

This has introduced breaking changes in 0.9 vs 0.8, and in speed up of
~100% ( 2.5 minutes instead of ~5minutes). The new reference data can be
found at
