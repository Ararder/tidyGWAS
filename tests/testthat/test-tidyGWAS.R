rs_merge_arch <- get_ref_data()
load(test_path("data/sumstats/test_sumstat.rds"))
load(test_path("data/sumstats/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))






# test initiate struct

test_that("Logging, name, and outdir owrks", {
  expect_no_error(test <- tidyGWAS(tbl = test_file, logfile = TRUE, name = "testrun", rs_merge_arch = rs_merge_arch))
  expect_no_error(test <- tidyGWAS(tbl = test_file, logfile = FALSE, name = "testrun", rs_merge_arch = rs_merge_arch))
  expect_no_error(test <- tidyGWAS(tbl = test_file, logfile = FALSE, rs_merge_arch = rs_merge_arch))
  expect_no_error(test <- tidyGWAS(tbl = test_file, logfile = TRUE, name = "auto-test", rs_merge_arch = rs_merge_arch))

})







# validate_snps----------------------------------------------------------

test_that("validate SNPs work", {

  struct <- initiate_struct(tbl = test_file, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(paste0(struct$filepaths$validate_snps))
  tmp <- validate_snps(struct, .filter_callback = .filter_callback)


})

test_that("validate_snps works, and detects failed parses of invalid RSID", {

  # add a row with a nonsensical RSID data
  tbl <- test_file |> dplyr::add_row(
    CHR = 1, POS = 101010,
    RSID = "XY_321332", EffectAllele = "G", OtherAllele = "T",
    EAF = 0.986, B =  -0.0262, SE = 0.0426, P = 0.539,
    CaseN = 106346, ControlN = 100000, INFO = 0.9)

  struct <- initiate_struct(tbl = tbl, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(paste0(struct$filepaths$validate_snps))
  expect_no_error(tmp <- validate_snps(struct, .filter_callback = .filter_callback))





})

test_that("testing validate stats", {


  struct <- initiate_struct(tbl = test_file, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(struct$filepaths$validate_stats)
  expect_no_error(tmp <- validate_stats(struct, .filter_callback = .filter_callback))


})



# full run ----------------------------------------------------------------



test_that("testing tidyGWAS without RSID", {
  skip("time consuming")
  bsgenome_objects <- get_bsgenome()
  tidyGWAS(
    tbl = dplyr::select(test_file, -RSID),
    bsgenome_objects = bsgenome_objects,
    rs_merge_arch = rs_merge_arch,
    logfile = TRUE,
    name = "full_test"
    )




})

test_that("testing tidyGWAS without CHR and POS", {
  skip("time consuming")

  tbl <- dplyr::select(test_file, -CHR, -POS) |>
    dplyr::tibble()
  tidyGWAS(
    tbl = tbl,
    bsgenome_objects = bsgenome_objects,
    rs_merge_arch = rs_merge_arch,
    logfile = TRUE,
    name = "full_test"
  )


})


test_that("Test time consuming functions", {
  skip("time consuming")

  tbl <-  dplyr::select(test_file, -RSID) |> dplyr::as_tibble()

  tidyGWAS(
    tbl = tbl,
    bsgenome_objects = bsgenome_objects,
    rs_merge_arch = rs_merge_arch,
    logfile = TRUE,
    name = "full_test"
  )



})






