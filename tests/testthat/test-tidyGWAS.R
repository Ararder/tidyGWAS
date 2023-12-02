







# validate_with_dbsnp -----------------------------------------------------

test_that("validate_with_dbsnp, all cols", {


  paths <- setup_pipeline_paths(tempfile())

  tbl <- dplyr::filter(flag_invalid_rsid(test_sumstat), !invalid_rsid)

  expect_no_error(validate_with_dbsnp(data_list, build = "NA", dbsnp_path = dbsnp_files,filepaths= paths))

  data_list$without_rsid <- NULL
  expect_no_error(validate_with_dbsnp(data_list, build = "NA", dbsnp_path = dbsnp_files, filepaths = paths))


})

test_that("validate_with_dbsnp, RSID", {

  paths <- setup_pipeline_paths(tempfile("tidyGWAS"))

  tmp <- dplyr::filter(flag_invalid_rsid(test_sumstat), !invalid_rsid) |>
    dplyr::select(-CHR, -POS)
  tmp_list <- list()
  tmp_list$main <- tmp


  expect_no_error(validate_with_dbsnp(tmp_list, build = "NA", dbsnp_path = dbsnp_files, paths))

})

test_that("validate_with_dbsnp, CHR and POS", {
  paths <- setup_pipeline_paths(tempfile("tidyGWAS"))
  tmp <- list()
  tmp$main <- dplyr::select(test_sumstat, -RSID) |>
    dplyr::mutate(CHR = as.character(CHR))


  expect_no_error(validate_with_dbsnp(tmp, build = "NA", dbsnp_path = dbsnp_files, filepaths = paths))

})





# tidyGWAS ----------------------------------------------------------------

test_that("can read in file from disk", {

  file <- withr::local_tempfile()
  arrow::write_csv_arrow(test_sumstat, file)

  expect_no_error(tmp <- tidyGWAS(tbl = file,  dbsnp_path = dbsnp_files))

})


test_that("Testing with CHR and POS", {


  expect_no_error(
    tidyGWAS(
      tbl = test_sumstat,
      dbsnp_path = dbsnp_files
    )
  )


  expect_no_error(
    tidyGWAS(
      tbl = dplyr::select(test_sumstat, -CHR, -POS),
      dbsnp_path = dbsnp_files
    )
  )

  expect_no_error(
    tidyGWAS(
      tbl = test_sumstat,
      dbsnp_path = dbsnp_files
    )
  )


})







test_that("Handles edge cases", {



    tfile <- flag_invalid_rsid(test_sumstat)
    tfile <- dplyr::mutate(tfile, CHR = dplyr::if_else(invalid_rsid, "50", CHR))
    tfile <- dplyr::mutate(tfile,  SE = dplyr::if_else(!invalid_rsid & CHR == "6", -50, SE))

    # edge case 1 - errors in both without_rsid and main
    expect_no_error(tidyGWAS(tfile, dbsnp_path = dbsnp_files, name = "edge-cases"))
    expect_no_error(tidyGWAS(pval_as_char_df, dbsnp_path = dbsnp_files, name = "edge-cases"))

    # test with indels
    expect_no_error(tidyGWAS(pval_as_char_df, dbsnp_path = dbsnp_files))
    expect_no_error(tidyGWAS(pval_as_char_df, dbsnp_path = dbsnp_files, keep_indels = FALSE))

    # handle where all rows are invalid_rsid
    tdf <- dplyr::filter(flag_invalid_rsid(test_sumstat), invalid_rsid)
    expect_no_error(tidyGWAS(tdf, dbsnp_path = dbsnp_files))

})



test_that("setup_pipeline_paths works", {

  expect_no_error(setup_pipeline_paths(tempfile()))

})


# -------------------------------------------------------------------------

test_that("write_finished_tidyGWAS works", {
  cleanup <- function(filepaths) {
    unlink(paste(filepaths$base, "tidyGWAS_hivestyle",sep="/"), recursive = TRUE)
    unlink(paste(filepaths$base, "cleaned_GRCh38.csv", sep="/"))
    unlink(paste(filepaths$base, "cleaned_GRCh38.csv.gz", sep="/"))
    unlink(paste(filepaths$base, "cleaned_GRCh38.parquet", sep="/"))

  }
  finished <- tidyGWAS(
    tbl = test_sumstat,
    logfile = TRUE,
    dbsnp_path = dbsnp_files,
  )
  filepaths <- setup_pipeline_paths(tempfile())
  cleanup(filepaths)
  expect_no_error(write_finished_tidyGWAS(finished, output_format = "hivestyle", outdir = tempdir(), filepaths = filepaths))
  cleanup(filepaths)
  expect_no_error(write_finished_tidyGWAS(finished, output_format = "csv", outdir = tempdir(), filepaths = filepaths))
  cleanup(filepaths)
  expect_no_error(write_finished_tidyGWAS(finished, output_format = "parquet", outdir = tempdir(), filepaths = filepaths))


})

