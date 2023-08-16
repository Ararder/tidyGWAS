test_that("arrow implementation works", {
  mock_arrow()

  chr_pos <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 37, by = "chr:pos")
  rsid <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_sumstat, 1:1000), build = 37, by = "rsid")

})





test_that("map_to_dbsnp_checks works", {
  expect_no_error(map_to_dbsnp_checks(test_sumstat, build = "37", by = "rsid", implementation = "bsgenome"))
  expect_no_error(map_to_dbsnp_checks(test_sumstat, build = "37", by = "rsid"))
  expect_no_error(map_to_dbsnp_checks(test_sumstat, build = "37", by = "chr:pos"))
  expect_no_error(map_to_dbsnp_checks(test_sumstat, build = "38", by = "chr:pos"))

})
