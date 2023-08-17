



# validate_sumstat----------------------------------------------------------

test_that("initiate_struct works", {
  filepaths <- setup_pipeline_paths("automated-testing", dbsnp_files)

  expect_no_error(struct <- initiate_struct(tbl = test_sumstat, filepaths = filepaths))

  # check that NAs are removed
  tmp <- dplyr::mutate(test_sumstat, B = dplyr::if_else(CHR == "21", NA_real_, B))
  n_na_rows <- nrow(dplyr::filter(test_sumstat,CHR == "21"))
  struct <- initiate_struct(tbl = tmp, filepaths = filepaths)
  expect_true(
    nrow(arrow::read_parquet(paste0(filepaths$removed_rows, "missing_values.parquet"))) == n_na_rows
  )
  # check that duplicates are removed
  tmp <- dplyr::bind_rows(test_file, dplyr::slice(test_file, c(1:10, 40:100, 10000:20000)))
  n_dups <- nrow(tmp) - nrow(test_sumstat)

  struct <- initiate_struct(tmp, filepaths = filepaths)
  expect_true(
    nrow(arrow::read_parquet(filepaths$duplicates)) == n_dups
  )

})

test_that("validate sumstat", {
  filepaths <- setup_pipeline_paths("automated-testing", dbsnp_files)

  struct <- initiate_struct(tbl = test_sumstat, filepaths)

  # if indels get passed initiate struct, should print weird alleles
  struct$sumstat <- pval_as_char_df
  struct$sumstat$rowid <- 1:nrow(pval_as_char_df)
  expect_message(tmp <- validate_sumstat(struct, verbose = FALSE, convert_p = 0))

})



test_that("validate SNPs even when 0 of invalid_rsids can be parsed", {

  tmp <- flag_invalid_rsid(test_sumstat) |>
    dplyr::mutate(RSID = dplyr::if_else(invalid_rsid, ".", RSID))

  struct <- initiate_struct(tbl = tmp, filepaths = setup_pipeline_paths("test", dbsnp_files))

  expect_no_error(tmp <- validate_sumstat(struct, convert_p = 0))


})

test_that("validate_snps works, and detects failed parses of invalid RSID", {

  # add a row with a nonsensical RSID data
  tbl <- test_sumstat |> dplyr::add_row(
    CHR = "1", POS = 101010,
    RSID = "XY_321332", EffectAllele = "G", OtherAllele = "T",
    EAF = 0.986, B =  -0.0262, SE = 0.0426, P = 0.539,
    CaseN = 106346, ControlN = 100000, INFO = 0.9)

  struct <- initiate_struct(tbl = tbl, filepaths = setup_pipeline_paths("test", dbsnp_files))

  expect_no_error(tmp <- validate_sumstat(struct, verbose=TRUE, convert_p = 2.225074e-308))


  # check that msg actually catches faulty rows
  struct$sumstat[10, "CHR"] <- "24"
  struct$sumstat[11, "POS"] <- -100L
  struct$sumstat[12, "EffectAllele"] <- "Y"
  struct$sumstat[13, "OtherAllele"] <- "X"

  expect_message(validate_sumstat(struct, convert_p = 2.225074e-308))




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
  expect_no_error(tmp <- validate_sumstat(struct, verbose=FALSE, convert_p = 2.225074e-308))


})





# validate_with_dbsnp -----------------------------------------------------

test_that("validate_with_dbsnp, all cols", {

  mock_arrow()

  struct <- initiate_struct(tbl = dplyr::filter(flag_invalid_rsid(test_sumstat), !invalid_rsid), filepaths = setup_pipeline_paths("test", dbsnp_files))

  expect_no_error(validate_with_dbsnp(struct, build = "NA"))


})

test_that("validate_with_dbsnp, RSID", {
  mock_arrow()

  tmp <- dplyr::select(test_sumstat, -CHR, -POS)

  struct <- initiate_struct(tbl = dplyr::filter(flag_invalid_rsid(test_sumstat), !invalid_rsid), filepaths = setup_pipeline_paths("test", dbsnp_files))
  expect_no_error(validate_with_dbsnp(struct, build = "NA"))

})

test_that("validate_with_dbsnp, CHR and POS", {
  mock_arrow()

  tmp <- dplyr::select(test_sumstat, -RSID) |>
    dplyr::mutate(CHR = as.character(CHR))

  struct <- initiate_struct(tbl = dplyr::filter(flag_invalid_rsid(test_sumstat), !invalid_rsid), filepaths = setup_pipeline_paths("test", dbsnp_files))
  expect_no_error(validate_with_dbsnp(struct, build = "NA"))

})





# tidyGWAS ----------------------------------------------------------------

test_that("can read in file from disk", {
  mock_arrow()
  file <- withr::local_tempfile()
  arrow::write_csv_arrow(test_file, file)

  tmp <- tidyGWAS(tbl = file, use_dbsnp = FALSE, dbsnp_path = dbsnp_files)

})


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


    tfile <- flag_invalid_rsid(test_sumstat)
    tfile <- dplyr::mutate(tfile, CHR = dplyr::if_else(invalid_rsid, "50", CHR))
    tfile <- dplyr::mutate(tfile,  SE = dplyr::if_else(!invalid_rsid & CHR == "6", -50, SE))

    # edge case 1 - errors in both without_rsid and main
    expect_no_error(tidyGWAS(tfile, dbsnp_path = dbsnp_files, name = "edge-cases"))

    # test with indels
    expect_no_error(tidyGWAS(pval_as_char_df, dbsnp_path = dbsnp_files))


})



test_that("setup_pipeline_paths works", {
  basename(withr::local_tempfile())
  expect_no_error(setup_pipeline_paths("testing", dbsnp_files))

})


# -------------------------------------------------------------------------

test_that("write_finished_tidyGWAS works", {
  cleanup <- function(filepaths) {
    unlink(paste(filepaths$base, "tidyGWAS_hivestyle",sep="/"), recursive = TRUE)
    unlink(paste(filepaths$base, "cleaned_GRCh38.csv", sep="/"))
    unlink(paste(filepaths$base, "cleaned_GRCh38.csv.gz", sep="/"))
    unlink(paste(filepaths$base, "cleaned_GRCh38.parquet", sep="/"))

  }
  mock_arrow()
  finished <- tidyGWAS(
    tbl = test_sumstat,
    logfile = TRUE,
    dbsnp_path = dbsnp_files,
    name = "test-write_finished_tidyGWAS"
  )
  filepaths <- setup_pipeline_paths("test-write_finished_tidyGWAS", dbsnp_files)
  cleanup(filepaths)
  expect_no_error(write_finished_tidyGWAS(finished, output_format = "hivestyle", outdir = tempdir(), filepaths = filepaths))
  cleanup(filepaths)
  expect_no_error(write_finished_tidyGWAS(finished, output_format = "csv", outdir = tempdir(), filepaths = filepaths))
  cleanup(filepaths)
  expect_no_error(write_finished_tidyGWAS(finished, output_format = "parquet", outdir = tempdir(), filepaths = filepaths))


})

