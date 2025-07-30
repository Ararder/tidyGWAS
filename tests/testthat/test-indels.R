test_that("multiplication works", {
  skip("long-range test")
  dbsnp_path <- "~/Downloads/dbSNP155"
  load(test_path("data/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))
  pval_as_char_df$rowid <- 1:nrow(pval_as_char_df)
  indel_strategy <- "qc"
  tbl <- pval_as_char_df
  filepaths <- setup_pipeline_paths(tempfile(), "raw")
  convert_p <- 0


  res <- detect_indels(
    tbl = tbl,
    indel_strategy = "qc",
    filepaths = filepaths,
    convert_p = convert_p
  )

})
