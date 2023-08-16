test_that("Repair_chr_pos works", {

  mock_arrow()
  test_sumstat$rowid <- 1:nrow(test_sumstat)
  tmp <- dplyr::select(test_sumstat, -CHR, -POS) |>
    flag_incorrect_rsid_format() |>
    dplyr::filter(!invalid_rsid)


  expect_no_error(repaired <- repair_chr_pos(tmp))

  # check that positions are imputed correctly
  in_dbsnp <- dplyr::filter(repaired, !no_dbsnp_entry) |>
    dplyr::select(CHR_37, POS_37, RSID, rowid)
  check <- dplyr::filter(test_sumstat, rowid %in% in_dbsnp$rowid) |>
    dplyr::select(CHR, POS, RSID, rowid)
  merged <- dplyr::inner_join(in_dbsnp, check, by ="rowid") |>
    dplyr::mutate(CHR = as.character(CHR))
  expect_equal(merged$POS_37, merged$POS)
  expect_equal(merged$CHR_37, merged$CHR)

})


test_that("Repair_rsid works", {
  mock_arrow()
  tmp <-
    dplyr::mutate(test_sumstat,  CHR = as.character(CHR), rowid = 1:nrow(test_sumstat)) |>
    dplyr::select(rowid, CHR, POS, RSID, EffectAllele, OtherAllele)


  expect_no_error(epaired <- repair_rsid(dplyr::select(tmp, -RSID)))

})



test_that("verify_chr_pos_rsid works", {
  mock_arrow()


  expect_no_error(verify_chr_pos_rsid(test_sumstat))
  expect_no_error(verify_chr_pos_rsid(test_sumstat, build = "37"))
  expect_error(verify_chr_pos_rsid(test_sumstat, build = 39))
  expect_no_error(verify_chr_pos_rsid(test_sumstat))


})






test_that("infer_build runs", {
  mock_arrow()
  expect_no_error(infer_build(test_sumstat))
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
