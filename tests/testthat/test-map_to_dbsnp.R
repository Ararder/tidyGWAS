test_that("arrow implementation works", {

  test_sumstat <- flag_invalid_rsid(test_sumstat) |> dplyr::filter(!invalid_rsid)
  expect_no_error(chr_pos <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 37, by = "chr:pos", dbsnp_path = dbsnp_files))
  expect_no_error(rsid <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 37, by = "rsid", dbsnp_path = dbsnp_files))

  expect_no_error(chr_pos <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 38, by = "chr:pos", dbsnp_path = dbsnp_files))
  expect_no_error(rsid <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 38, by = "rsid", dbsnp_path = dbsnp_files))

})


test_that("map_to_dbsnp works", {
  expect_no_error(map_to_dbsnp(test_sumstat, build = "37", by = "chr:pos", dbsnp_path = dbsnp_files))

  expect_error(map_to_dbsnp(dplyr::tibble(.rows = 1), build = "37", by = "chr:pos", dbsnp_path = dbsnp_files))
})
