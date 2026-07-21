test_that("extra input columns can be retained", {
  input <- test_sumstat |>
    dplyr::slice_head(n = 100) |>
    dplyr::mutate(
      source_label = paste0("source-", rowid),
      source_integer = as.integer(rowid * 2),
      source_all_na = NA_character_
    )
  output_dir <- file.path(withr::local_tempdir(), "tidyGWAS-output")

  result <- tidyGWAS(
    tbl = input,
    dbsnp_path = dbsnp_path,
    build = "37",
    repair_cols = FALSE,
    keep_columns = TRUE,
    output_format = "parquet",
    output_dir = output_dir
  )

  expected_rows <- match(result$rowid, input$rowid)
  expect_identical(result$source_label, input$source_label[expected_rows])
  expect_identical(result$source_integer, input$source_integer[expected_rows])
  expect_identical(result$source_all_na, input$source_all_na[expected_rows])

  saved <- arrow::read_parquet(
    file.path(output_dir, "tidyGWAS_cleaned.parquet")
  )
  expect_true(
    all(c("source_label", "source_integer", "source_all_na") %in% colnames(saved))
  )
  expect_identical(saved$source_label, result$source_label)
})


test_that("extra input columns are dropped by default", {
  input <- test_sumstat |>
    dplyr::slice_head(n = 100) |>
    dplyr::mutate(source_label = paste0("source-", rowid))

  result <- tidyGWAS(
    tbl = input,
    dbsnp_path = dbsnp_path,
    build = "37",
    repair_cols = FALSE,
    output_dir = tempfile()
  )

  expect_false("source_label" %in% colnames(result))
})


test_that("extra columns cannot overwrite columns created by tidyGWAS", {
  main <- dplyr::tibble(rowid = 1:2, indel = c(FALSE, TRUE))
  extra_tbl <- dplyr::tibble(rowid = 1:2, indel = c("input", "input"))

  expect_error(
    restore_extra_columns(main, extra_tbl, "indel"),
    "tidyGWAS creates output column"
  )
})


test_that("rowid must be present, complete, and unique", {
  expect_error(
    validate_keep_columns(dplyr::tibble(value = 1:2)),
    "must contain a.*rowid"
  )
  expect_error(
    validate_keep_columns(dplyr::tibble(rowid = c(1L, NA_integer_))),
    "cannot contain missing values"
  )
  expect_error(
    validate_keep_columns(dplyr::tibble(rowid = c(1L, 1L))),
    "must contain unique values"
  )

  valid <- dplyr::tibble(rowid = 1:2, value = c("a", "b"))
  expect_identical(validate_keep_columns(valid), valid)
})


test_that("reserved tidyGWAS output names are rejected before processing", {
  input <- dplyr::tibble(rowid = 1:2)

  expect_error(
    validate_keep_columns(input, c("indel", "discrep_freq")),
    "reserves these output column names"
  )
  expect_identical(validate_keep_columns(input, "cohort"), input)
})
