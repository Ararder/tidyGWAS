test_that("dups are removed", {
  filepaths <- setup_pipeline_paths(tempfile())
  fp <- tempfile()
  expect_no_error(remove_duplicates(tbl = test_sumstat, filepath = ))
  expect_true(!file.exists(fp))

  # add some duplications
  tmp <- dplyr::bind_rows(test_sumstat[1:10,], test_sumstat)
  tmp$rowid <- 1:nrow(tmp)
  remove_duplicates(tbl = tmp, filepath = fp)
  expect_true(file.exists(fp))
  expect_equal(nrow(arrow::read_parquet(fp)), 10)

})




test_that("NA rows are removed", {

  filepaths <- setup_pipeline_paths(tempfile())

  test_sumstat$RSID <- dplyr::if_else(test_sumstat$CHR == "6", NA_character_, test_sumstat$RSID)
  expect_no_error(remove_rows_with_na(test_sumstat, columns = "CHR", filepath = filepaths$removed_missing_critical))


})



test_that("indels are removed", {

  filepaths <- setup_pipeline_paths(tempfile())


  expect_no_error(detect_indels(pval_as_char_df, FALSE, filepaths, convert_p = 0))
  expect_no_error(detect_indels(pval_as_char_df, TRUE, filepaths, convert_p = 0))

  #



})


