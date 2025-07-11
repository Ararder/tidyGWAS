


# tidyGWAS ----------------------------------------------------------------



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

test_that("POS column is removed if exist in indels",{

    res <- tidyGWAS(
      tbl = pval_as_char_df,
      dbsnp_path = dbsnp_path
    )
    expect_false("POS" %in% colnames(res))

})


test_that("Handles only RSID and no non-rsid values", {
  tmp_test <- test_sumstat
  tmp_test$CHR <- NULL
  tmp_test$POS <- NULL
  tmp_test1 <- dplyr::filter(tmp_test, stringr::str_detect(RSID, "rs"))


  expect_no_error(
    tidyGWAS(
      tbl = tmp_test1,
      dbsnp_path = dbsnp_path
    )
  )

  expect_no_error(
    tidyGWAS(
      tbl = tmp_test,
      dbsnp_path = dbsnp_path
    )
  )



})

test_that("Handles duplications when chr-pos in RSID maps to existing row", {
  tmp_test <- test_sumstat

  add <- dplyr::slice(tmp_test,1:15) |>
    dplyr::mutate(RSID = paste0(CHR, ":", POS)) |>
    dplyr::select(-rowid)

  add2 <- dplyr::slice(tmp_test,1:15) |>
    dplyr::mutate(RSID = paste0("27", ":", POS)) |>
    dplyr::select(-rowid)


  tmp_test2 <- dplyr::bind_rows(add2, tmp_test) |> dplyr::select(-rowid, -CHR,-POS)
  tmp_test <- dplyr::bind_rows(add, tmp_test) |> dplyr::select(-rowid)
  tmp_test$CHR <- NULL
  tmp_test$POS <- NULL



  expect_no_error(
    tidyGWAS(
      tbl = tmp_test,
      dbsnp_path = dbsnp_path
    )
  )

  expect_no_error(
    tidyGWAS(
      tbl = tmp_test2,
      dbsnp_path = dbsnp_path
    )
  )



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



test_that("add EAF flag", {

  test_tbl <- dplyr::mutate(test_sumstat, EAF = EAF + 0.2)
  main <- tidyGWAS(
    tbl = test_tbl,
    dbsnp_path = dbsnp_path,
    flag_discrep_freq = "AMR"
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

test_that("repair allele frequency works", {
  skip("requires local")
  dbsnp_path <- "~/Downloads/dbSNP155/"

  expect_no_error(tidyGWAS(
    dplyr::select(tbl, -EAF),
    dbsnp_path = dbsnp_path,
    impute_freq = "EUR"
  ))


  expect_no_error(tidyGWAS(
    dplyr::select(tbl, -EAF),
    dbsnp_path = dbsnp_path,
    impute_freq = "AFR"
  ))

})



test_that("allele freq flag check works", {


  tidyGWAS(
    tbl,
    dbsnp_path = dbsnp_path,
    flag_discrep_freq = "AMR"
  )

})



test_that("Filter EAF works", {


  tidyGWAS(
    tbl = dplyr::select(test_sumstat, -EAF),
    dbsnp_path = dbsnp_path
    # min_EAF = 0.25
  )

})




