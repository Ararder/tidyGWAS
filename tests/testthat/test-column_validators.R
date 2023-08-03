

## validate RSID -----------------------------------------------------------

test_that("validate RSIDs does not error", {
  # create rowid
  tmp <- dplyr::mutate(test_file, rowid = 1:nrow(test_file))

  only_correct_rsid <- flag_incorrect_rsid_format(tmp) |> dplyr::filter(!invalid_rsid)
  only_incorrect_rsid <- flag_incorrect_rsid_format(tmp) |> dplyr::filter(invalid_rsid)
  expect_message(validate_rsid(tmp), regexp = "Succeded in parsing 307 out of 307 rows*")

  expect_no_error(validate_rsid(only_correct_rsid))


})





## validate B ------------------------------------------------------------


test_that("Validate_b can detect a mislabelled odds-ratio column", {
  expect_message(
    validate_b(dplyr::mutate(test_file, B = exp(B))),
    regexp = "*indicating that B has been mislabelled*"
  )

  expect_no_error(validate_b(test_file))


})


## validate P -------------------------------------------------------------------------


test_that("Validate P catches pvalues that are converted to character", {
  expect_message(validate_P(pval_as_char_df))

})




## validate SE -------------------------------------------------------------------------

test_that("Validate SE work", {
  expect_message(validate_se(pval_as_char_df))

  pval_as_char_df$SE[3] <- -1
  expect_true(nrow(validate_se(pval_as_char_df) |> dplyr::filter(invalid_SE)) == 1)

})


## validate N -------------------------------------------------------------------------

test_that("Validate N runs without error", {
  expect_message(validate_n(pval_as_char_df))

})





## validate EA and OA ------------------------------------------------------


test_that("validate_effect_cols detects and prints alleles outside ACGT", {



  expect_message(validate_ea_oa(pval_as_char_df))


})


# validate POS

test_that("Validate POS runs", {

  expect_message(validate_pos(pval_as_char_df))


})


# validate CHR ------------------------------------------------------------

test_that("Validat CHR runs, and fixes UCSC format", {


  before <- test_file$CHR
  tmp <- dplyr::mutate(test_file, CHR = dplyr::if_else(EffectAllele == "T", as.character(CHR), paste0("chr", CHR)))
  after <- validate_chr(tmp)
  expect_equal(as.character(before), after$CHR)

})



# validate Z ------------------------------------------------------------

test_that("Z", {
  tmp <- dplyr::mutate(test_file, Z = B/SE)
  # unreasonable large Z
  tmp[100, ]$Z <- 500

  # should detect extremely large Z
  expect_message(validate_z(tmp), regexp =  "WARNING*")
  # otherwise should be fine
  expect_no_error(validate_z(tmp[-100, ]))


})

