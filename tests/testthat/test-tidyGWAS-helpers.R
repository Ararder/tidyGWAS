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


test_that("select_correct_columns work", {
  tmp <- dplyr::select(test_sumstat, CHR, EffectAllele, OtherAllele)
  expect_error(
    select_correct_columns(
      tbl = tmp
      )
    )


})



test_that("correct error printing for column_names", {
  tmp <- dplyr::rename(test_sumstat, CHROM = CHR, A1 = EffectAllele, A2 = OtherAllele)
  names <- list(
    CHR = "CHROM",
    EffectAllele = "A1",
    OtherAllele = "A2",
    INFO = "pop",
    CaseN = "cCAse"
  )

  expect_error(tidyGWAS(
    tmp,
    dbsnp_path,
    column_names = names
    )
  )

})

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

