test_that("Repair_chr_pos works", {

  mock_arrow()

  tmp <- dplyr::select(test_sumstat, -CHR, -POS) |>
    flag_invalid_rsid() |>
    dplyr::filter(!invalid_rsid)


  expect_no_error(repaired <- repair_chr_pos(tmp, dbsnp_path = dbsnp_files))


})


test_that("Repair_rsid works", {
  mock_arrow()
  tmp <-
    dplyr::mutate(test_sumstat,  CHR = as.character(CHR), rowid = 1:nrow(test_sumstat)) |>
    dplyr::select(rowid, CHR, POS, RSID, EffectAllele, OtherAllele)


  expect_no_error(repair_rsid(dplyr::select(tmp, -RSID), dbsnp_path = dbsnp_files))

})



test_that("verify_chr_pos_rsid works", {
  mock_arrow()


  expect_no_error(verify_chr_pos_rsid(test_sumstat, dbsnp_path = dbsnp_files))

  expect_error(verify_chr_pos_rsid(test_sumstat, build = 39))



})






test_that("infer_build runs", {
  mock_arrow()
  expect_no_error(infer_build(test_sumstat, dbsnp_path = dbsnp_files))
})

test_that("flatten_dbsnp runs", {
  expect_no_error(flatten_dbsnp(b38))
})


test_that("check_incompat_alleles runs", {
 tmp <- dplyr::mutate(test_sumstat, rowid = 1:nrow(test_sumstat))

 expect_no_error(check_incompat_alleles(tmp, flatten_dbsnp(b38)))
})

test_that("make callback runs", {
  expect_no_error(make_callback(withr::local_tempfile()))
})
