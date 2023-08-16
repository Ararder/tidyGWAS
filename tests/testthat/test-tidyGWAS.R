



# validate_sumstat----------------------------------------------------------

test_that("initiate_struct works", {
  filepaths <- setup_pipeline_paths("automated-testing", dbsnp_files)
  expect_no_error(struct <- initiate_struct(tbl = test_sumstat, filepaths = filepaths))

})

test_that("validate sumstat", {
  filepaths <- setup_pipeline_paths("automated-testing", dbsnp_files)

  struct <- initiate_struct(tbl = test_sumstat, filepaths)

  # if indels get passed initiate struct, should print weird alleles
  struct$sumstat <- pval_as_char_df
  struct$sumstat$rowid <- 1:nrow(pval_as_char_df)
  expect_message(tmp <- validate_sumstat(struct, verbose = FALSE))

})



test_that("validate SNPs even when 0 of invalid_rsids can be parsed", {

  tmp <- flag_incorrect_rsid_format(test_sumstat) |>
    dplyr::mutate(RSID = dplyr::if_else(invalid_rsid, ".", RSID))

  struct <- initiate_struct(tbl = tmp, filepaths = setup_pipeline_paths("test", dbsnp_files))

  expect_no_error(tmp <- validate_sumstat(struct))


})

test_that("validate_snps works, and detects failed parses of invalid RSID", {

  # add a row with a nonsensical RSID data
  tbl <- test_sumstat |> dplyr::add_row(
    CHR = "1", POS = 101010,
    RSID = "XY_321332", EffectAllele = "G", OtherAllele = "T",
    EAF = 0.986, B =  -0.0262, SE = 0.0426, P = 0.539,
    CaseN = 106346, ControlN = 100000, INFO = 0.9)

  struct <- initiate_struct(tbl = tbl, filepaths = setup_pipeline_paths("test", dbsnp_files))

  expect_no_error(tmp <- validate_sumstat(struct, verbose=TRUE))


  # check that msg actually catches faulty rows
  struct$sumstat[10, "CHR"] <- "24"
  struct$sumstat[11, "POS"] <- -100L
  struct$sumstat[12, "EffectAllele"] <- "Y"
  struct$sumstat[13, "OtherAllele"] <- "X"

  expect_message(validate_sumstat(struct))




})

test_that("test that validate_sumstat catches errors in columns", {

  # setup
  struct <- initiate_struct(tbl = test_sumstat,filepaths = setup_pipeline_paths("test", dbsnp_files))



  struct$sumstat[100, "B"] <- Inf
  struct$sumstat[101, "P"] <- -4
  struct$sumstat[102, "SE"] <- 0
  struct$sumstat[100, "SE"] <- 0
  struct$sumstat[103, "N"] <- 0
  struct$sumstat[104, "EAF"] <- 1
  # five less rows should exist after validation
  expect_no_error(tmp <- validate_sumstat(struct, verbose=FALSE))


})





# validate_with_dbsnp -----------------------------------------------------

test_that("validate_with_dbsnp, all cols", {

  mock_arrow()
  struct <- initiate_struct(tbl = test_sumstat, filepaths = setup_pipeline_paths("test", dbsnp_files))
  expect_no_error(validate_with_dbsnp(struct))


})

test_that("validate_with_dbsnp, RSID", {
  mock_arrow()

  tmp <- dplyr::select(test_sumstat, -CHR, -POS)

  struct <- initiate_struct(tbl = tmp, filepaths = setup_pipeline_paths("test", dbsnp_files))
  expect_no_error(validate_with_dbsnp(struct))

})

test_that("validate_with_dbsnp, CHR and POS", {
  mock_arrow()

  tmp <- dplyr::select(test_sumstat, -RSID) |>
    dplyr::mutate(CHR = as.character(CHR))

  struct <- initiate_struct(tbl = tmp, filepaths = setup_pipeline_paths("test", dbsnp_files))
  expect_no_error(validate_with_dbsnp(struct))

})





# tidyGWAS ----------------------------------------------------------------









test_that("Testing with CHR and POS", {
  mock_arrow()


  expect_no_error(
    tidyGWAS(
      tbl = dplyr::select(test_sumstat, -RSID),
      logfile = TRUE,
      dbsnp_path = dbsnp_files
    )
  )


  expect_no_error(
    tidyGWAS(
      tbl = dplyr::select(test_sumstat, -CHR, -POS),
      logfile = FALSE,
      dbsnp_path = dbsnp_files
    )
  )

  expect_no_error(
    tidyGWAS(
      tbl = test_sumstat,
      logfile = FALSE,
      dbsnp_path = dbsnp_files
    )
  )

  expect_no_error(tidyGWAS(test_sumstat, dbsnp_path = dbsnp_files))

})







test_that("Handles edge cases", {

    mock_arrow()


    tfile <- flag_incorrect_rsid_format(test_sumstat)
    tfile <- dplyr::mutate(tfile, CHR = dplyr::if_else(invalid_rsid, "50", CHR))
    tfile <- dplyr::mutate(tfile,  SE = dplyr::if_else(!invalid_rsid & CHR == "6", -50, SE))

    # edge case 1 - errors in both without_rsid and main
    expect_no_error(tidyGWAS(tfile, dbsnp_path = dbsnp_files))

    # test with indels
    expect_no_error(tidyGWAS(pval_as_char_df, dbsnp_path = dbsnp_files))


})



test_that("setup_pipeline_paths works", {
  basename(withr::local_tempfile())
  expect_no_error(setup_pipeline_paths("testing", dbsnp_files))

})
