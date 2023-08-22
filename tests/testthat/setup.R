load(test_path("data/test_sumstat.rds"))
load(test_path("data/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))

pval_as_char_df$rowid <- 1:nrow(pval_as_char_df)

test_sumstat <- dplyr::tibble(test_sumstat)
test_sumstat$rowid <- 1:nrow(test_sumstat)
test_sumstat$CHR <- as.character(test_sumstat$CHR)
dbsnp_files <- test_path("fixtures/dbSNP155")
dbsnp_path <- test_path("fixtures/dbSNP155")


# default args
name <- "testing"
convert_p <- 0
verbose <- FALSE
keep_indels <- TRUE




