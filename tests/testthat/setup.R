load(test_path("data/test_sumstat.rds"))
test_sumstat <- dplyr::tibble(test_sumstat)
test_sumstat$CHR <- as.character(test_sumstat$CHR)
load(test_path("data/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))
pval_as_char_df$rowid <- 1:nrow(pval_as_char_df)
dbsnp_path <- test_path("fixtures/dbSNP155")



data_list <- list()
data_list$main <- dplyr::filter(flag_invalid_rsid(test_sumstat), !invalid_rsid)
data_list$without_rsid <- dplyr::filter(flag_invalid_rsid(test_sumstat), invalid_rsid) |>
  dplyr::select(-RSID)

# default args
logfile <- FALSE
name <- "testing"
convert_p <- 0
verbose <- FALSE
keep_indels <- TRUE
overwrite = FALSE
build <- "NA"
tbl <- test_sumstat
indel_strategy <- "keep"
# outdir = paste0(tempdir(), "/",stringr::str_replace_all(date(), pattern = c(" "="_", ":"="_")))
# output_dir = paste0(tempdir(), "/",stringr::str_replace_all(date(), pattern = c(" "="_", ":"="_")))



# reference:
#   - title: Munging summary statistics made easy
#
# - subtitle: dbSNP155
# desc: >
#   Functions that interact with the dbSNP155 database.
# contents:
#   - tidyGWASls
# - infer_build
# - update_ids
#
# - subtitle: Cleaning
# desc: >
#   Functions that detect and clean up common issues in GWAS summary statistics.
# contents:
#   - flag_indels
# - detect_indels
# - flag_duplicates
# - flag_invalid_rsid
# - remove_duplicates
# - remove_rows_with_na
# - select_correct_columns
#
# - subtitle: Validation
# desc: >
#   Functions that validate GWAS columns
# contents:
#   - validate_columns
# - validate_rsid
# - validate_sumstat

