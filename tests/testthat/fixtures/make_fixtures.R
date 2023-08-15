# rs_merge_arch <- tidyGWAS::rs_merge_arch
# save(rs_merge_arch, file = test_path("fixtures/rs_merge_arch.rds"))
#
# b38 <- map_to_dbsnp(dplyr::tibble(test_file), by = "rsid", build = 38)
# b37 <- map_to_dbsnp(dplyr::tibble(test_file), by = "rsid", build = 37)
#
#
# save(b38, file = test_path("fixtures/b38.rds"))
# save(b37, file = test_path("fixtures/b37.rds"))

# setup local rs_merge_arch ------------------------------------------------


# load(test_path("data/sumstats/test_sumstat.rds"))
# load(test_path("data/sumstats/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))
#
# temp <- dplyr::bind_rows(
#   dplyr::select(test_file, RSID),
#   dplyr::select(pval_as_char_df, RSID)
# ) |>
#   dplyr::mutate(rowid = 1:127715) |>
#   flag_rsid_history()
#
#
# dplyr::filter(temp, !is.na(old_RSID)) |>
#   dplyr::select(old_RSID, RSID) |>
#   arrow::write_parquet(test_path("fixtures/RsMergeArch.parquet"))



# setup_local_dbsnp -------------------------------------------------------
local_dbsnp <- function() {
  first <- dplyr::select(test_file, CHR, POS, RSID)
  second <- dplyr::select(pval_as_char_df, CHR, POS, RSID)
  second$CHR <- as.character(second$CHR)

  data <- dplyr::bind_rows(first, second) |>
    dplyr::distinct() |>
    dplyr::tibble()

  grch37_rsid <- map_to_dbsnp(tbl = data, build = 37, by = "rsid")
  grch37_chr_pos <- map_to_dbsnp(tbl = data, build = 37, by = "chr:pos")

  final <- dplyr::bind_rows(grch37_chr_pos,grch37_rsid) |>
    dplyr::distinct() |>
    dplyr::mutate(RSID = as.integer(stringr::str_sub(RSID, start = 3)))


  final$CHR <- as.character(final$CHR)
  arrow::write_dataset(dplyr::group_by(final, CHR), test_path("fixtures/dbSNP/GRCh37"), compression = "gzip")


  grch38_rsid <- map_to_dbsnp(tbl = data, build = 38, by = "rsid")
  grch38_chr_pos <- map_to_dbsnp(tbl = data, build = 38, by = "chr:pos")


  final <- dplyr::bind_rows(grch38_rsid,grch38_chr_pos) |>
    dplyr::distinct() |>
    dplyr::mutate(RSID = as.integer(stringr::str_sub(RSID, start = 3)))
  final$CHR <- as.character(final$CHR)
  arrow::write_dataset(dplyr::group_by(final, CHR), test_path("fixtures/dbSNP/GRCh38"), compression = "gzip")


}
