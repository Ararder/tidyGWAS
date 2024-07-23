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
outdir = paste0(tempdir(), "/",stringr::str_replace_all(date(), pattern = c(" "="_", ":"="_")))
output_dir = paste0(tempdir(), "/",stringr::str_replace_all(date(), pattern = c(" "="_", ":"="_")))



