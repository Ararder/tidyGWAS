#

test_that("map_to_dbsnp works", {
  expect_no_error(map_to_dbsnp(test_sumstat, build = "37", dbsnp_path = dbsnp_path))

  expect_error(map_to_dbsnp(dplyr::tibble(.rows = 1), build = "37", dbsnp_path = dbsnp_path))
})

test_that("flag no_dbsnp_entry works", {
  

  ll <- repair_ids(tbl, repair = "rsid", dbsnp_path = dbsnp_path)
  expect_equal(sum(ll$no_dbsnp_entry), 21)


})

test_that("Repair_ids works", {



  tmp <- dplyr::select(test_sumstat, -CHR, -POS) |>
    flag_invalid_rsid() |>
    dplyr::filter(!invalid_rsid)


  expect_no_error(repaired <- repair_ids(tmp, repair = "pos", dbsnp_path = dbsnp_path))

  tmp <-
    dplyr::mutate(test_sumstat,  CHR = as.character(CHR), rowid = 1:nrow(test_sumstat)) |>
    dplyr::select(rowid, CHR, POS, RSID, EffectAllele, OtherAllele)


  expect_no_error(
    repair_ids(
      dplyr::select(tmp, -RSID),
      repair = "rsid",
      dbsnp_path = dbsnp_path)
    )

})









test_that("infer_build runs", {

  expect_no_error(infer_build(test_sumstat, dbsnp_path = dbsnp_path))
})


test_that("make callback runs", {
  expect_no_error(make_callback(withr::local_tempfile()))
})
