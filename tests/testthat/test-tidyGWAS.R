



# validate_snps----------------------------------------------------------


test_that("validate SNPs work", {

  struct <- initiate_struct(tbl = test_file, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(paste0(struct$filepaths$validate_snps))
  expect_no_error(
    tmp <- validate_snps(struct, .filter_callback = .filter_callback, verbose = FALSE)
    )

  # if indels get passed initiate struct, should print weird alleles
  struct$sumstat <- pval_as_char_df
  struct$sumstat$rowid <- 1:nrow(pval_as_char_df)
  expect_message(tmp <- validate_snps(struct, .filter_callback = .filter_callback, verbose = FALSE))

})



test_that("validate SNPs even when 0 of invalid_rsids can be parsed", {
  tmp <- flag_incorrect_rsid_format(test_file) |>
    dplyr::mutate(RSID = dplyr::if_else(invalid_rsid, ".", RSID))

  struct <- initiate_struct(tbl = tmp, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(paste0(struct$filepaths$validate_snps))
  expect_no_error(tmp <- validate_snps(struct, .filter_callback = .filter_callback,verbose=FALSE))


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
  expect_no_error(tmp <- validate_snps(struct, .filter_callback = .filter_callback, verbose=TRUE))


  # check that msg actually catches faulty rows
  struct$sumstat[10, "CHR"] <- "24"
  struct$sumstat[11, "POS"] <- -100L
  struct$sumstat[12, "EffectAllele"] <- "Y"
  struct$sumstat[13, "OtherAllele"] <- "X"





})

test_that("testing validate stats", {

  # setup
  struct <- initiate_struct(tbl = test_file, rs_merge_arch = rs_merge_arch, filepaths = setup_pipeline_paths("test"))
  .filter_callback = make_callback(struct$filepaths$validate_stats)


  expect_no_error(
    tmp <- validate_stats(struct, .filter_callback = .filter_callback, verbose=TRUE)
  )





  struct$stats[100, "B"] <- Inf
  struct$stats[101, "P"] <- -4
  struct$stats[102, "SE"] <- 0
  struct$stats[100, "SE"] <- 0
  struct$stats[103, "N"] <- 0
  struct$stats[104, "EAF"] <- 1
  expect_no_error(
    tmp <- validate_stats(struct, .filter_callback = .filter_callback, verbose=FALSE)
    )


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
  mock_dbsnp()
  bsgenome <- bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")

  expect_no_error(
      tidyGWAS(
      tbl = test_file,
      rs_merge_arch = rs_merge_arch,
      bsgenome_objects = bsgenome,
      name = "full",
      verbose = FALSE
  ))





})


test_that("testing without CHR and POS", {
  skip("Full length test to run on longleaf")
  # mock_dbsnp()
  # bsgenome <- bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")

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

test_that("Testing without CHR and POS", {
  mock_dbsnp()
  bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")


  expect_no_error(

    tidyGWAS(
      tbl = dplyr::select(test_file, -CHR, -POS),
      bsgenome_objects = bsgenome_objects,
      rs_merge_arch = rs_merge_arch,
      logfile = FALSE,
      name = "chr_pos"
    )

  )

})



test_that("Testing with minimum input of parameters", {
  mock_dbsnp()
  bsgenome_objects <- list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")


  expect_no_error(tidyGWAS(test_file))

})







