# Package index

## tidyGWAS

### Cleaning summary statistics

- [`tidyGWAS()`](https://ararder.github.io/tidyGWAS/reference/tidyGWAS.md)
  : Execute validation and quality control of GWAS summmary statistics

### Interact with GWAS catalog and Open Targets

- [`from_gwas_catalog()`](https://ararder.github.io/tidyGWAS/reference/from_gwas_catalog.md)
  : Download summary statistics from GWAS catalog
- [`from_gwas_catalog_region()`](https://ararder.github.io/tidyGWAS/reference/from_gwas_catalog_region.md)
  : Query a specific region of interest for a using a gwas catalog
  study_id
- [`get_open_targets_cs()`](https://ararder.github.io/tidyGWAS/reference/get_open_targets_cs.md)
  : Query Open targets for all credible sets containing the variant
- [`check_rest_avail()`](https://ararder.github.io/tidyGWAS/reference/check_rest_avail.md)
  : Check if API access is available for a GWAS catalog study

### Meta-analysis

- [`meta_analyse()`](https://ararder.github.io/tidyGWAS/reference/meta_analyse.md)
  : Improved meta-analysis using tidyGWAS:ed files

- [`create_lake()`](https://ararder.github.io/tidyGWAS/reference/create_lake.md)
  : Create a data lake in hivestyle format

- [`deprec_meta_analyze()`](https://ararder.github.io/tidyGWAS/reference/deprec_meta_analyze.md)
  :

  Perform meta-analysis of GWAS summary statistics datasets cleaned by
  tidyGWAS This function is depreciated. Use
  [`meta_analyse()`](https://ararder.github.io/tidyGWAS/reference/meta_analyse.md)

- [`align_to_ref()`](https://ararder.github.io/tidyGWAS/reference/align_to_ref.md)
  : Align EffectAllele to always be the reference genome allele

- [`by_chrom()`](https://ararder.github.io/tidyGWAS/reference/by_chrom.md)
  : Meta analyse a chromosome, one at a time!

### Ancestry estimation

- [`ancestry_comp()`](https://ararder.github.io/tidyGWAS/reference/ancestry_comp.md)
  : Estimate ancestry composition from allele frequencies

### Helpful functions

- [`flag_duplicates()`](https://ararder.github.io/tidyGWAS/reference/flag_duplicates.md)
  : Find all rows which are part of a set of duplicated rows
- [`flag_indels()`](https://ararder.github.io/tidyGWAS/reference/flag_indels.md)
  : Detect "indels" in GWAS summary statistics
- [`flag_invalid_rsid()`](https://ararder.github.io/tidyGWAS/reference/flag_invalid_rsid.md)
  : Detect entries that are not valid rsID's in GWAS summary statistics
- [`strand_flip()`](https://ararder.github.io/tidyGWAS/reference/strand_flip.md)
  : Strand flip alleles
- [`validate_rsid()`](https://ararder.github.io/tidyGWAS/reference/validate_rsid.md)
  : Validate format of the RSID column in a GWAS summary statistics file
- [`infer_build()`](https://ararder.github.io/tidyGWAS/reference/infer_build.md)
  : Infer what genome build a GWAS summary statistics file is on.
- [`repair_ids()`](https://ararder.github.io/tidyGWAS/reference/repair_ids.md)
  : Augment a data.frame with information from dbSNP
- [`repair_stats()`](https://ararder.github.io/tidyGWAS/reference/repair_stats.md)
  : Repair statistics column in a GWAS summary statistics tibble
- [`validate_sumstat()`](https://ararder.github.io/tidyGWAS/reference/validate_sumstat.md)
  : Validate statistics columns in a GWAS summary statistics file
