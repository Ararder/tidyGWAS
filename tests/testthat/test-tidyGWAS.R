


# tidyGWAS ----------------------------------------------------------------

test_that("Can update column names", {
  column_map <- list(
    CHR = "chrom",
    POS = "position",
    RSID = "SNP",
    EffectAllele = "A1",
    OtherAllele = "A2"
  )
  tbl <- dplyr::rename(
    tbl,
    chrom = CHR,
    position = POS,
    SNP = RSID,
    A1 = EffectAllele,
    A2 = OtherAllele
  )

  expect_no_error(
    update_column_names(tbl, column_map, CaseN = NULL, ControlN = NULL, N = NULL)
  )


})


test_that("Can parse file from disk with different column names", {
  column_map <- list(
    CHR = "chrom",
    POS = "position",
    RSID = "SNP",
    EffectAllele = "A1",
    OtherAllele = "A2"
  )
  tbl <- dplyr::rename(
    tbl,
    chrom = CHR,
    position = POS,
    SNP = RSID,
    A1 = EffectAllele,
    A2 = OtherAllele
  )

  file <- withr::local_tempfile()
  arrow::write_csv_arrow(tbl, file)

  expect_no_error(
    tidyGWAS(
      tbl = file,
      dbsnp_path = dbsnp_path,
      column_names = column_map,
      output_dir = tempfile()
    )
  )


})

test_that("can read in file from disk", {

  file <- withr::local_tempfile()
  arrow::write_csv_arrow(test_sumstat, file)

  expect_no_error(tmp <- tidyGWAS(tbl = file,  dbsnp_path = dbsnp_path, output_dir = tempfile()))

})

test_that("indels handled correctly", {
  cleaned <- tidyGWAS(
    tbl = pval_as_char_df,
    dbsnp_path,
    indel_strategy = "remove"
  )
  expect_true(nrow(cleaned) < nrow(pval_as_char_df) -1500)

  cleaned <- tidyGWAS(
    tbl = pval_as_char_df,
    dbsnp_path,
    indel_strategy = "keep"
  )
  expect_true(nrow(cleaned) == nrow(pval_as_char_df) -161)

})
test_that("Testing with CHR and POS", {


  expect_no_error(
    tidyGWAS(
      tbl = test_sumstat,
      dbsnp_path = dbsnp_path
    )
  )


  expect_no_error(
    tidyGWAS(
      tbl = dplyr::select(test_sumstat, -CHR, -POS),
      dbsnp_path = dbsnp_path
    )
  )


  expect_no_error(
    tidyGWAS(
      tbl = dplyr::select(test_sumstat, -RSID),
      dbsnp_path = dbsnp_path
    )
  )



})


test_that("regression test", {
    cleaned <- tidyGWAS(
      tbl = test_sumstat,
      dbsnp_path = dbsnp_path
    )
  expect_true(
    nrow(cleaned) == 99979
  )
})


test_that("minimum number of columns", {
  expect_no_error(
    cleaned <- tidyGWAS(
      tbl = dplyr::select(test_sumstat, RSID, EffectAllele, OtherAllele),
      dbsnp_path = dbsnp_path
    )
  )
  expect_no_error(
    cleaned <- tidyGWAS(
      tbl = dplyr::select(test_sumstat, CHR,POS, EffectAllele, OtherAllele),
      dbsnp_path = dbsnp_path
    )
  )

  expect_error(
    cleaned <- tidyGWAS(
      tbl = dplyr::select(test_sumstat, CHR,EffectAllele, OtherAllele),
      dbsnp_path = dbsnp_path
    )
  )
})




test_that("Handles edge cases", {
  skip()


  tfile <- flag_invalid_rsid(test_sumstat)
  tfile <- dplyr::mutate(tfile, CHR = dplyr::if_else(invalid_rsid, "50", CHR))
  tfile <- dplyr::mutate(tfile,  SE = dplyr::if_else(!invalid_rsid & CHR == "6", -50, SE))

  # edge case 1 - errors in both without_rsid and main
  expect_no_error(tidyGWAS(tfile, dbsnp_path = dbsnp_path))
  expect_no_error(tidyGWAS(pval_as_char_df, dbsnp_path))

  # test with indels
  expect_no_error(tidyGWAS(pval_as_char_df, dbsnp_path = dbsnp_path))
  expect_no_error(tidyGWAS(pval_as_char_df, dbsnp_path = dbsnp_path, indel_strategy = "remove"))

  # handle where all rows are invalid_rsid
  tdf <- dplyr::filter(flag_invalid_rsid(test_sumstat), invalid_rsid)
  expect_no_error(tidyGWAS(tdf, dbsnp_path = dbsnp_path))

})



test_that("setup_pipeline_paths works", {

  expect_no_error(setup_pipeline_paths(tempfile()))

})


# -------------------------------------------------------------------------

test_that("write_finished_tidyGWAS works", {

  finished <- tidyGWAS(
    tbl = test_sumstat,
    logfile = TRUE,
    dbsnp_path = dbsnp_path,
    output_dir = tempfile()
  )
  filepaths <- setup_pipeline_paths(tempfile("testing"))

  unlink(filepaths$base, recursive = TRUE)
  expect_no_error(write_finished_tidyGWAS(finished, output_format = "hivestyle",  filepaths = filepaths))


})

