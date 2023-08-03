

# initiate struct ---------------------------------------------------------



test_that("Logging, name, and outdir owrks", {
  expect_no_error(test <- tidyGWAS(tbl = test_file, logfile = TRUE, name = "testrun", rs_merge_arch = rs_merge_arch))
  expect_no_error(test <- tidyGWAS(tbl = test_file, logfile = FALSE, name = "testrun", rs_merge_arch = rs_merge_arch))
  expect_no_error(test <- tidyGWAS(tbl = test_file, logfile = FALSE, rs_merge_arch = rs_merge_arch))

})




# validate_snps----------------------------------------------------------


test_that("validate SNPs work", {

  struct <- initiate_struct(tbl = test_file, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(paste0(struct$filepaths$validate_snps))
  expect_no_error(tmp <- validate_snps(struct, .filter_callback = .filter_callback))


})


test_that("validate SNPs even when 0 of invalid_rsids can be parsed", {
  tmp <- flag_incorrect_rsid_format(test_file) |>
    dplyr::mutate(RSID = dplyr::if_else(invalid_rsid, ".", RSID))

  struct <- initiate_struct(tbl = tmp, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(paste0(struct$filepaths$validate_snps))
  expect_no_error(tmp <- validate_snps(struct, .filter_callback = .filter_callback))


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
  bsgenome <- bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")

  struct <- initiate_struct(tbl = test_file, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  expect_no_error(validate_with_dbsnp(struct, bsgenome_objects = bsgenome, make_callback(struct$filepaths$validate_with_dbsnp)))

})

test_that("validate_with_dbsnp, RSID", {
  mock_dbsnp()
  bsgenome <- bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")
  tmp <- dplyr::select(test_file, -CHR, -POS)

  struct <- initiate_struct(tbl = tmp, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  expect_no_error(validate_with_dbsnp(struct, bsgenome_objects = bsgenome, .filter_callback = make_callback(struct$filepaths$validate_with_dbsnp)))

})

test_that("validate_with_dbsnp, CHR and POS", {
  mock_dbsnp()
  bsgenome <- bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")
  tmp <- dplyr::select(test_file, -RSID) |>
    dplyr::mutate(CHR = as.character(CHR))

  struct <- initiate_struct(tbl = tmp, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  expect_no_error(validate_with_dbsnp(struct, bsgenome_objects = bsgenome, make_callback(struct$filepaths$validate_with_dbsnp)))

})


# tidyGWAS ----------------------------------------------------------------


test_that("testing tidyGWAS with RSID, CHR and POS", {
  skip("covr github actions fail")
  mock_dbsnp()
  bsgenome <- bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")
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
  bsgenome <- bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")

  tmp_file <- data.table::fread("~/shared/gwas_sumstats/sumstats/default/default/raw/ukb_phase1to3_rest_edge_full_dec21_2019_pheno87.fastGWA.gz") |>
    dplyr::rename(RSID = SNP, EffectAllele = A1, OtherAllele = A2, EAF = AF1, B = BETA) |>
    dplyr::tibble()
    # dplyr::slice_sample(n = 1000000)

  tmp_file <- data.table::fread("~/shared/gwas_sumstats/sumstats/insomnia/insomnia/raw/insomnia.cc.tsv.gz") |>
    dplyr::rename(RSID = SNP, POS = BP, EffectAllele = A1, OtherAllele = A2, EAF = FREQ, B = BETA) |>
    dplyr::tibble()

  lk <- tidyGWAS(
    tbl = tmp_file,
    bsgenome_objects = bs,
    rs_merge_arch = rs_merge_arch,
    name = "no_chr_pos"
  )


})


test_that("Testing without RSID", {
  mock_dbsnp()
  bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")


  expect_no_error(

  tidyGWAS(
    tbl = dplyr::select(test_file, -RSID),
    bsgenome_objects = bsgenome_objects,
    rs_merge_arch = rs_merge_arch,
    logfile = TRUE,
    name = "no_rsid"
  )

  )

})







