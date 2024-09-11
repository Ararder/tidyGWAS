

test_that("meta-analysis works with variable number of columns", {
  ea <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, 10)
  dset <- setup_test_data_meta_analysis(EAF = c(0.1, 0.3, 0.3), N = n)
  
  res <- meta_analyze(dset)

  expect_true(!"INFO" %in% colnames(res))
  
  # EA should correspond to weighted by sample size
  expected <- sum(ea * n) / sum(n)
  expect_true(res$EAF == expected)


})



test_that("meta-analysis handles missing values in calculation of INFO/EAF", {
  ea <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, NA)
  info <- c(0.6, 0.85, 1.3)
  dset <- setup_test_data_meta_analysis(EAF = ea, N = n, INFO = info)
  
  res <- meta_analyze(dset)
  
  # EA should correspond to weighted by sample size
  expected <- sum(ea * n, na.rm = T) / sum(n, na.rm =T)
  expect_true(res$EAF == expected)

  # now for INFO
  expected_i <- sum(info * n, na.rm = T) / sum(n, na.rm = T)
  expect_true(res$INFO == expected_i)


})




test_that("meta-analysis works with EAF and not INFO", {
  vals <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, 10)
  
  dset <- setup_test_data_meta_analysis(EAF = vals, N = n)
  res <- meta_analyze(dset)
  
  expect_true(!"INFO" %in% colnames(res))
  expected <- sum(vals * n) / sum(n)
  expect_true(res$EAF == expected)



})



test_that("meta-analysis works with INFO and not EAF", {
  vals <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, 10)
  
  dset <- setup_test_data_meta_analysis(INFO = vals, N = n)
  res <- meta_analyze(dset)
  
  expect_true(!"EAF" %in% colnames(res))
  expected <- sum(vals * n) / sum(n)
  expect_true(res$INFO == expected)


})


test_that("meta-analysis works with variable number of columns", {
  
  dset <- setup_test_data_meta_analysis()
  expect_no_error(res <- meta_analyze(dset))
  


})
