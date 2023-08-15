test_that("arrow implementation works", {
  skip("time consuming")

  chr_pos <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_file, 1:1000), build = 37, by = "chr:pos")

  tictoc::tic("arrow with rsid")
  rsid <- map_to_dbsnp_arrow(tbl = dplyr::slice(test_file, 1:1000), build = 37, by = "rsid")
  tictoc::toc()
})



test_that("bsgenome implementation works", {
  skip("time consuming")
  old_test <- dplyr::filter(old_test, CHR != "2")
  tictoc::tic("bsgenome")
  rsid <- map_to_dbsnp_bsgenome(tbl = test_file, build = 37, by = "rsid")
  tictoc::toc()
  chr_pso <- map_to_dbsnp_bsgenome(tbl = test_file, build = 37, by = "chr:pos")

})

test_that("map_to_dbsnp_checks works", {
  expect_no_error(map_to_dbsnp_checks(test_file, build = "37", by = "rsid", implementation = "bsgenome"))
  expect_no_error(map_to_dbsnp_checks(test_file, build = "37", by = "rsid"))
  expect_no_error(map_to_dbsnp_checks(test_file, build = "37", by = "chr:pos"))
  expect_no_error(map_to_dbsnp_checks(test_file, build = "38", by = "chr:pos"))

})
