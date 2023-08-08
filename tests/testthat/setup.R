load(test_path("fixtures/rs_merge_arch.rds"))
load(test_path("data/sumstats/test_sumstat.rds"))
load(test_path("data/sumstats/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))
load(test_path("fixtures/b38.rds"))
load(test_path("fixtures/b37.rds"))
test_file$P <- as.character(test_file$P)
add_indels <- dplyr::filter(flag_indels(pval_as_char_df), indel) |> dplyr::mutate(CaseN = 50, ControlN = 30, INFO = 0.8) |>
  dplyr::select(-N, -indel)
test_file <- dplyr::bind_rows(test_file, add_indels)


test_file$CHR <- as.character(test_file$CHR)
test_file <- dplyr::tibble(test_file)
test_file$rowid <- 1:nrow(test_file)


