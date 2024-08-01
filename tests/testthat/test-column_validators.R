

## validate RSID -----------------------------------------------------------

test_that("validate RSIDs does not error", {
  # what happens if validate_rid cannot parse format?
  filepaths <- setup_pipeline_paths(tempfile())





  # some invalid rsids, no CHR or POS

  expect_no_error(validate_rsid(dplyr::select(test_sumstat, -CHR, -POS), filepaths$removed_invalid_rsid))

  # some invalid rsids, no CHR or POS, some failed
  tmp <- dplyr::select(test_sumstat, -CHR, -POS)
  tmp <- dplyr::mutate(tmp, RSID = dplyr::if_else(B > 0, ".", RSID))
  expect_no_error(validate_rsid(tmp, filepaths$removed_invalid_rsid))

  # all invalid
  tmp <- dplyr::select(test_sumstat, -CHR, -POS)
  tmp <- dplyr::mutate(tmp, RSID = ".")
  expect_no_error(validate_rsid(tmp, filepaths$removed_invalid_rsid))




})

test_that("validate RSIDs does not error", {
  # create rowid
  tmp <- dplyr::mutate(test_sumstat, rowid = 1:nrow(test_sumstat)) |>
    dplyr::select(-CHR, -POS)

  only_correct_rsid <- flag_invalid_rsid(tmp) |> dplyr::filter(!invalid_rsid)
  only_incorrect_rsid <- flag_invalid_rsid(tmp) |> dplyr::filter(invalid_rsid)
  expect_no_error(validate_rsid(tmp))

  expect_no_error(validate_rsid(only_correct_rsid))
  expect_no_error(validate_rsid(tmp))



})

test_that("Validate_columns works", {


  expect_no_error(validate_sumstat(tbl, convert_p = 0))

})




test_that("Validate_columns works", {

  check <- colnames(test_sumstat)[colnames(test_sumstat) %in% stats_cols]
  check <- check[!check %in% c("CaseN", "ControlN", "INFO")]
  test_sumstat <- dplyr::mutate(test_sumstat, P = dplyr::if_else(CHR == 6, -3, P))

  expect_no_error(
    for(c in check) tmp <- validate_columns(test_sumstat,col = c, convert_p = 0)
  )



})

test_that("validate sumstat", {
  filepaths <- setup_pipeline_paths(tempfile())

  pval_as_char_df$rowid <- 1:nrow(pval_as_char_df)

  expect_message(tmp <- validate_sumstat(pval_as_char_df, convert_p = 0))

})







test_that("validate_columns catches issues", {
  expect_message(
    validate_columns(dplyr::mutate(test_sumstat, B = exp(B)), "B"),
    regexp = "*indicating that B has been mislabelled*"
  )

})







## validate EA and OA ------------------------------------------------------

test_that("validate_effect_cols detects and prints alleles outside ACGT", {


  expect_message(validate_columns(pval_as_char_df, "EffectAllele"))
  expect_message(validate_columns(pval_as_char_df, "OtherAllele"))

})








# validate CHR ------------------------------------------------------------

test_that("Validate CHR runs, and fixes UCSC format", {

  before <- test_sumstat$CHR
  tmp <- dplyr::mutate(test_sumstat, CHR = dplyr::if_else(EffectAllele == "T", as.character(CHR), paste0("chr", CHR)))
  after <- validate_columns(tmp, "CHR")
  expect_equal(as.character(before), after$CHR)

})



# validate Z ------------------------------------------------------------

test_that("Z", {
  tmp <- dplyr::mutate(test_sumstat, Z = B/SE)
  # unreasonable large Z
  tmp[100, ]$Z <- 500

  # should detect extremely large Z
  expect_message(validate_columns(tmp, "Z"), regexp =  "WARNING*")
  # otherwise should be fine
  expect_no_error(validate_columns(tmp[-100, ], "Z"))


})



# Can handle missing values in all columns

test_that("validate_columns can handle missing values", {

  tmp <- test_sumstat |>
    dplyr::mutate(
      Z = B/SE,
      dplyr::across(c("B","SE","EAF", "INFO", "P", "CaseN", "ControlN", "Z"), \(x) dplyr::if_else(CHR == 6, NA_real_, x))
    )

  expect_no_error(validate_columns(tmp, "B"))
  expect_no_error(validate_columns(tmp, "SE"))
  expect_no_error(validate_columns(tmp, "EAF"))
  expect_no_error(validate_columns(tmp, "P", convert_p=0))
  expect_no_error(validate_columns(tmp, "Z"))



})

