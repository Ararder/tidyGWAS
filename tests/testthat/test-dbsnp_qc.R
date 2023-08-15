test_that("Repair_chr_pos works", {

  mock_arrow()
  test_file$rowid <- 1:nrow(test_file)
  tmp <- dplyr::select(test_file, -CHR, -POS) |>
    flag_incorrect_rsid_format() |>
    dplyr::filter(!invalid_rsid)


  expect_no_error(repaired <- repair_chr_pos(tmp))

  # check that positions are imputed correctly
  in_dbsnp <- dplyr::filter(repaired, !no_dbsnp_entry) |>
    dplyr::select(CHR_37, POS_37, RSID, rowid)
  check <- dplyr::filter(test_file, rowid %in% in_dbsnp$rowid) |>
    dplyr::select(CHR, POS, RSID, rowid)
  merged <- dplyr::inner_join(in_dbsnp, check, by ="rowid") |>
    dplyr::mutate(CHR = as.character(CHR))
  expect_equal(merged$POS_37, merged$POS)
  expect_equal(merged$CHR_37, merged$CHR)

})


test_that("Repair_rsid works", {
  mock_arrow()
  tmp <-
    dplyr::mutate(test_file,  CHR = as.character(CHR), rowid = 1:nrow(test_file)) |>
    dplyr::select(rowid, CHR, POS, RSID, EffectAllele, OtherAllele)


  expect_no_error(epaired <- repair_rsid(dplyr::select(tmp, -RSID)))

})



test_that("verify_chr_pos_rsid works", {
  mock_arrow()
  tmp <-
    dplyr::mutate(test_file,  CHR = as.character(CHR), rowid = 1:nrow(test_file)) |>
    dplyr::select(rowid, CHR, POS, RSID, EffectAllele, OtherAllele)

  expect_no_error(verify_chr_pos_rsid(tmp, build = "38"))
  expect_no_error(verify_chr_pos_rsid(tmp, build = "37"))
  expect_error(verify_chr_pos_rsid(tmp, build = 39))
  expect_no_error(verify_chr_pos_rsid(tmp))


})







test_that("infer_build runs", {
  mock_arrow()
  expect_no_error(infer_build(test_file))
})

test_that("flatten_dbsnp runs", {
  expect_no_error(flatten_dbsnp(b38))
})


test_that("qc_with_dbsnp runs", {
 tmp <- dplyr::mutate(test_file, rowid = 1:nrow(test_file))

 expect_no_error(check_incompat_alleles(tmp, flatten_dbsnp(b38)))
})

test_that("make callback runs", {
  expect_no_error(make_callback(withr::local_tempfile()))
})
