load(test_path("data/test_sumstat.rds"))
test_sumstat <- dplyr::tibble(test_sumstat)
test_sumstat$CHR <- as.character(test_sumstat$CHR)
load(test_path("data/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))
pval_as_char_df$rowid <- 1:nrow(pval_as_char_df)
dbsnp_path <- test_path("fixtures/dbSNP155")



# default args
logfile <- FALSE
name <- "testing"
convert_p <- 0
verbose <- FALSE
overwrite = FALSE
build <- "NA"
tbl <- test_sumstat
output_dir = tempfile()
# tmp <- readr::read_tsv("~/Desktop/local.tsv")
#
#
# column_names = list(
#   CHR = "#CHROM",
#   RSID = "ID",
#   EffectAllele = "A1",
#   OtherAllele = "REF",
#   EAF = "A1_FREQ",
#   SE = "LOG(OR)_SE",
#   Z = "Z_STAT",
#   N = "OBS_CT"
# )
#
# tmp
# #
# # res <- tidyGWAS::tidyGWAS(
# #   tmp,
# #   dbsnp_path = "~/Downloads/dbSNP155",
# #   column_names = list(
# #     CHR = "#CHROM",
# #     RSID = "ID",
# #     EffectAllele = "A1",
# #     OtherAllele = "REF",
# #     EAF = "A1_FREQ",
# #     SE = "LOG(OR)_SE",
# #     Z = "Z_STAT",
# #     N = "OBS_CT"
# #   ))
