load(test_path("fixtures/rs_merge_arch.rds"))
load(test_path("data/sumstats/test_sumstat.rds"))
load(test_path("data/sumstats/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))
test_file <- dplyr::tibble(test_file)
test_file$CHR <- as.character(test_file$CHR)






# initiate struct ---------------------------------------------------------



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


test_that("validate SNPs even when 0 of invalid_rsids can be parsed", {
  tmp <- flag_incorrect_rsid_format(test_file) |>
    dplyr::mutate(RSID = dplyr::if_else(invalid_rsid, ".", RSID))

  struct <- initiate_struct(tbl = tmp, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(paste0(struct$filepaths$validate_snps))
  tmp <- validate_snps(struct, .filter_callback = .filter_callback)


})

test_that("validate_snps works, and detects failed parses of invalid RSID", {

  # add a row with a nonsensical RSID data
  tbl <- test_file |> dplyr::add_row(
    CHR = "1", POS = 101010,
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

# validate_with_dbsnp -----------------------------------------------------

test_that("validate_with_dbsnp, all cols", {
  mock_dbsnp()
  bsgenome <- get_bsgenome()

  struct <- initiate_struct(tbl = test_file, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  expect_no_error(validate_with_dbsnp(struct, bsgenome_objects = bsgenome, make_callback(struct$filepaths$validate_with_dbsnp)))

})

test_that("validate_with_dbsnp, RSID", {
  mock_dbsnp()
  bsgenome <- get_bsgenome()
  tmp <- dplyr::select(test_file, -CHR, -POS)

  struct <- initiate_struct(tbl = tmp, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  expect_no_error(validate_with_dbsnp(struct, bsgenome_objects = bsgenome, .filter_callback = make_callback(struct$filepaths$validate_with_dbsnp)))

})

test_that("validate_with_dbsnp, CHR and POS", {
  mock_dbsnp()
  bsgenome <- get_bsgenome()
  tmp <- dplyr::select(test_file, -RSID) |>
    dplyr::mutate(CHR = as.character(CHR))

  struct <- initiate_struct(tbl = tmp, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  expect_no_error(validate_with_dbsnp(struct, bsgenome_objects = bsgenome, make_callback(struct$filepaths$validate_with_dbsnp)))

})


# tidyGWAS ----------------------------------------------------------------



test_that("testing tidyGWAS with RSID, CHR and POS", {
  skip("covr github actions fail")
  mock_dbsnp()
  bsgenome <- get_bsgenome()
  expect_no_error(
    tidyGWAS(
    tbl = test_file,
    rs_merge_arch = rs_merge_arch,
    bsgenome_objects = bsgenome,
    name = "full"
    ))




})

test_that("testing without CHR and POS", {
  skip("covr github actions fail")
  mock_dbsnp()
  bsgenome <- get_bsgenome()

  tbl <- dplyr::select(test_file, -CHR, -POS) |>
    dplyr::tibble()
  expect_no_error(

  tidyGWAS(
    tbl = tbl,
    bsgenome_objects = bsgenome_objects,
    rs_merge_arch = rs_merge_arch,
    name = "no_chr_pos"
  )

  )

})


test_that("Testing without RSID", {
  mock_dbsnp()
  bsgenome <- get_bsgenome()

  tbl <- dplyr::select(test_file, -RSID) |>
    dplyr::tibble()
  expect_no_error(

  tidyGWAS(
    tbl = tbl,
    bsgenome_objects = bsgenome_objects,
    rs_merge_arch = rs_merge_arch,
    logfile = FALSE,
    name = "no_rsid"
  )

  )

})






