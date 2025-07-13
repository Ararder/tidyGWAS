
test_that("meta-analysis works with variable number of columns", {


  ea <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, 10)
  dset <- setup_test_data_meta_analysis(EAF = c(0.1, 0.3, 0.3), N = n)

  res <- meta_analyze(dset)
  res2 <- meta_analyse2(dset, ref = "REF_37")
  expect_no_error(res2 <- meta_analyse2(dset, ref = "REF_37"))

  expect_true(!"INFO" %in% colnames(res))

  # EA should correspond to weighted by sample size
  expected <- sum(ea * n) / sum(n)
  expect_true(res$EAF == expected)


})



test_that("regression test", {
  ea <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, NA)
  info <- c(0.6, 0.85, 1.3)
  dset <- setup_test_data_meta_analysis(EAF = ea, N = n, INFO = info)


  res2 <- meta_analyse2(dset)
  res <- meta_analyze(dset)

  expect_true(res$B == res2$B)
  expect_true(res$P == res2$P)


})



test_that("meta-analysis handles missing values in calculation of INFO/EAF", {
  ea <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, NA)
  info <- c(0.6, 0.85, 1.3)
  dset <- setup_test_data_meta_analysis(EAF = ea, N = n, INFO = info)


  res2 <- meta_analyse2(dset)
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
  res <- meta_analyze2(dset)

  expect_true(!"EAF" %in% colnames(res))
  expected <- sum(vals * n) / sum(n)
  expect_true(res$INFO == expected)


})


test_that("meta-analysis works with variable number of columns", {

  dset <- setup_test_data_meta_analysis()
  expect_no_error(res <- meta_analyze(dset))



})



test_that("benchmark new version", {
  skip()
  two <- arrow::open_dataset("~/sumstats/iq/tidyGWAS_hivestyle")
  one <- arrow::open_dataset("~/sumstats/aers_plus/tidyGWAS_hivestyle", schema = arrow::schema(two))
  dset <- c(one,two)



  by = c("RSID", "EffectAllele", "OtherAllele")
  ref = "REF_37"

  tictoc::tic("old")
  k <- meta_analyze_by_chrom(dset, chrom = "2", by = c("POS_38", "RSID", "EffectAllele", "OtherAllele"), ref =  "REF_37")
  tictoc::toc()

  tictoc::tic("old")
  k2 <- meta_analyze(dset, ref="REF_37", impl = "old")
  tictoc::toc()

  tictoc::tic("new")
  k2 <- meta_analyze(dset, ref="REF_37", impl = "new")
  tictoc::toc()




})


# -------------------------------------------------------------------------
# complex tests
ds <- setup_test_data_meta_analysis(
  test_sumstat,
  n_variants_per_chr = 5,
  n_traits = 10,
  frac_missing_N = 0.1,
  frac_missing_EAF = 0.1,
  frac_missing_INFO = 0.1,
  drop_N_for_traits = 1,
  drop_EAF_for_traits = FALSE,
  drop_INFO_for_traits = FALSE,
  outdir = tempfile()
)


test_that("meta-analysis correctly calculates N", {
  df <- ds |> dplyr::collect() |>

  eaf_truth <-
    df |>
    dplyr::filter(RSID == "rs80311739") |>
    dplyr::select(RSID, EAF, EffectiveN) |>
    dplyr::mutate(
      t = EAF * EffectiveN,
      n_w_eaf = dplyr::if_else(is.na(EAF), 0, EffectiveN)
      ) |>
    dplyr::summarise(
      EAF = sum(t, na.rm=TRUE) / sum(n_w_eaf,na.rm=TRUE)
      )

  res <- meta_analyse2(ds)

  dplyr::filter(res, RSID == "rs80311739")




})

