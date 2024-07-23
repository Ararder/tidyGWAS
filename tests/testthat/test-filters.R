# parse_tbl
test_that("validate parse_tbl", {
  tmp_file <- withr::local_tempfile(fileext = ".csv")
  test_sumstat <- dplyr::tibble(test_sumstat)
  tbl <- arrow::write_csv_arrow(test_sumstat, tmp_file)

  # can read csv file
  expect_equal(test_sumstat, tbl)

  expect_equal(tbl, parse_tbl(test_sumstat)$tbl)

  # can handle data.frame
  expect_equal(tbl, parse_tbl(as.data.frame(test_sumstat))$tbl)

  # errors if not passed a filepath or data.frame
  expect_error(parse_tbl("list"))

  #
  expect_error(parse_tbl())


})

test_that("select_correct_columns work", {
  munged <- select_correct_columns(test_sumstat)

  # should not remove any columns from test_file
  expect_true(all(colnames(test_sumstat) %in% colnames(munged)))

  # if all columns are removed, should throw an error
  expect_error(
    purrr::set_names(test_sumstat, c("a","b","c","d","e","f","g","h", "x", "y", "z","p")) |>
      select_correct_columns()
  )

  # can handle that study_n is not passed
  # one column with all NAs should not result in all rows being removed
  tmp <- dplyr::mutate(test_sumstat, EAF = NA_real_)
  expect_true(nrow(select_correct_columns(tmp)) == nrow(test_sumstat))

})



test_that("dups are removed", {
  filepaths <- setup_pipeline_paths(tempfile())

  expect_no_error(remove_duplicates(tbl = test_sumstat, filepaths = filepaths))
  expect_true(!file.exists(filepaths$removed_duplicates))

  # add some duplications
  tmp <- dplyr::bind_rows(test_sumstat[1:10,], test_sumstat)
  tmp$rowid <- 1:nrow(tmp)
  remove_duplicates(tbl = tmp, filepaths = filepaths)
  expect_true(file.exists(filepaths$removed_duplicates))
  expect_equal(nrow(arrow::read_parquet(filepaths$removed_duplicates)), 10)

})




test_that("NA rows are removed", {

  filepaths <- setup_pipeline_paths(tempfile())

  test_sumstat$RSID <- dplyr::if_else(test_sumstat$CHR == "6", NA_character_, test_sumstat$RSID)
  expect_no_error(remove_rows_with_na(test_sumstat, filepaths))


})



test_that("indels are removed", {

  filepaths <- setup_pipeline_paths(tempfile())


  expect_no_error(detect_indels(pval_as_char_df, FALSE, filepaths, convert_p = 0))
  expect_no_error(detect_indels(pval_as_char_df, TRUE, filepaths, convert_p = 0))

  #



})


