# tidyGWAS 0.9.0

Major revision from earlier versions. The reference data format has been updated:

1.  Compression switched from gzip to snappy (better speed, and less installation issues)
2.  GRCh37 and GRCh38 have been merged into a single file, requiring only one merge to get CHR, POS and RSID for both builds.

This has introduced breaking changes in 0.9 vs 0.8, and in speed up of \~100% ( 2.5 minutes instead of \~5minutes). The new reference data can be found at
