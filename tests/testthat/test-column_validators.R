# source(test_path("setup.R"))

## validate RSID -----------------------------------------------------------

test_that("validate RSIDs does not error", {
  # create rowid
  tmp <- dplyr::mutate(test_file, rowid = 1:nrow(test_file))

  only_correct_rsid <- flag_incorrect_rsid_format(tmp) |> dplyr::filter(!invalid_rsid)
  only_incorrect_rsid <- flag_incorrect_rsid_format(tmp) |> dplyr::filter(invalid_rsid)
  expect_no_error(validate_rsid(tmp))

  expect_no_error(validate_rsid(only_correct_rsid))

  # test without chr and pos
  validate_rsid(dplyr::select(test_file, -CHR, -POS))


})



test_that("Validate_columns works", {

  check <- colnames(test_file)[colnames(test_file) %in% stats_cols]
  check <- check[!check %in% c("CaseN", "ControlN", "INFO")]
  test_file <- dplyr::mutate(test_file, P = dplyr::if_else(CHR == 6, "-3", P))

  expect_no_error(
    for(c in check) tmp <- validate_columns(test_file,col = c, verbose = FALSE)
  )


})






test_that("validate_columns catches issues", {
  expect_message(
    validate_columns(dplyr::mutate(test_file, B = exp(B)), "B"),
    regexp = "*indicating that B has been mislabelled*"
  )

  expect_message(validate_columns(pval_as_char_df, "P"))

})







## validate EA and OA ------------------------------------------------------


test_that("validate_effect_cols detects and prints alleles outside ACGT", {

  expect_message(validate_columns(pval_as_char_df, "EffectAllele"))
  expect_message(validate_columns(pval_as_char_df, "OtherAllele"))

})





# validate CHR ------------------------------------------------------------

test_that("Validat CHR runs, and fixes UCSC format", {


  before <- test_file$CHR

  tmp <- dplyr::mutate(test_file, CHR = dplyr::if_else(EffectAllele == "T", as.character(CHR), paste0("chr", CHR)))
  after <- validate_columns(tmp, "CHR")
  before <- dplyr::if_else(before == "23", "X", before)
  expect_equal(as.character(before), after$CHR)

})



# validate Z ------------------------------------------------------------

test_that("Z", {
  tmp <- dplyr::mutate(test_file, Z = B/SE)
  # unreasonable large Z
  tmp[100, ]$Z <- 500

  # should detect extremely large Z
  expect_message(validate_columns(tmp, "Z"), regexp =  "WARNING*")
  # otherwise should be fine
  expect_no_error(validate_columns(tmp[-100, ], "Z"))


})

