# test_that("arrow implementation works", {
#
#   test_sumstat <- flag_invalid_rsid(test_sumstat) |> dplyr::filter(!invalid_rsid)
#   expect_no_error(chr_pos <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 37, by = "chr:pos", dbsnp_path = dbsnp_files))
#   expect_no_error(rsid <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 37, by = "rsid", dbsnp_path = dbsnp_files))
#
#   expect_no_error(chr_pos <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 38, by = "chr:pos", dbsnp_path = dbsnp_files))
#   expect_no_error(rsid <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 38, by = "rsid", dbsnp_path = dbsnp_files))
#
# })
#

test_that("map_to_dbsnp works", {
  expect_no_error(map_to_dbsnp(test_sumstat, build = "37", by = "chr:pos", dbsnp_path = dbsnp_files))

  expect_error(map_to_dbsnp(dplyr::tibble(.rows = 1), build = "37", by = "chr:pos", dbsnp_path = dbsnp_files))
})

test_that("Repair_chr_pos works", {



  tmp <- dplyr::select(test_sumstat, -CHR, -POS) |>
    flag_invalid_rsid() |>
    dplyr::filter(!invalid_rsid)


  expect_no_error(repaired <- repair_chr_pos(tmp, dbsnp_path = dbsnp_files))


})


test_that("Repair_rsid works", {

  tmp <-
    dplyr::mutate(test_sumstat,  CHR = as.character(CHR), rowid = 1:nrow(test_sumstat)) |>
    dplyr::select(rowid, CHR, POS, RSID, EffectAllele, OtherAllele)


  expect_no_error(repair_rsid(dplyr::select(tmp, -RSID), dbsnp_path = dbsnp_files))

})








test_that("infer_build runs", {

  expect_no_error(infer_build(test_sumstat, dbsnp_path = dbsnp_files))
})



test_that("make callback runs", {
  expect_no_error(make_callback(withr::local_tempfile()))
})
