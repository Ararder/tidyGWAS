

test_that("meta-analysis works with variable number of columns", {


  ea <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, 10)
  dset <- setup_test_data_meta_analysis(EAF = c(0.1, 0.3, 0.3), N = n)

  res <- meta_analyse(dset)

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


  res <- meta_analyse(dset)



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
  res <- meta_analyse(dset)

  expect_true(!"INFO" %in% colnames(res))
  expected <- sum(vals * n) / sum(n)
  expect_true(res$EAF == expected)



})



test_that("meta-analysis works with INFO and not EAF", {
  vals <- c(0.1, 0.3, 0.3)
  n <- c(10, 15, 10)

  dset <- setup_test_data_meta_analysis(INFO = vals, N = n)
  res <- meta_analyse(dset)

  expect_true(!"EAF" %in% colnames(res))
  expected <- sum(vals * n) / sum(n)
  expect_true(res$INFO == expected)


})






# -------------------------------------------------------------------------
# complex tests
ds <- setup_test_data_meta_analysis2(
  test_sumstat,
  n_variants_per_chr = 5,
  n_traits = 10,
  frac_missing_N = 0.1,
  frac_missing_EAF = 0.1,
  frac_missing_INFO = 0.1,
  drop_N_for_traits = FALSE,
  drop_EAF_for_traits = FALSE,
  drop_INFO_for_traits = FALSE,
  outdir = tempfile()
)


test_that("meta-analysis correctly handles missing vallues in INFO and EAF", {
  df <- ds |> dplyr::collect()

  eaf_truth <-
    df |>
      dplyr::group_by(RSID, EffectAllele,OtherAllele) |>
      dplyr::select(RSID, EffectAllele, OtherAllele,EAF, N, INFO) |>
      dplyr::mutate(
        t = EAF * N,
        i = INFO * N,
        n_w_eaf = dplyr::if_else(is.na(EAF), 0, N),
        n_w_info = dplyr::if_else(is.na(INFO), 0, N)
        ) |>
      dplyr::summarise(
        EAF = sum(t, na.rm=TRUE) / sum(n_w_eaf,na.rm=TRUE),
        INFO = sum(i, na.rm = TRUE) / sum(n_w_info, na.rm = TRUE),
        N = sum(N, na.rm = TRUE)
        ) |>
      dplyr::arrange(RSID, EffectAllele, OtherAllele)

  res <- meta_analyse(ds) |> dplyr::arrange(RSID, EffectAllele, OtherAllele)

  expect_equal(res$EAF, eaf_truth$EAF, tolerance = 1e-5)
  expect_equal(res$INFO, eaf_truth$INFO, tolerance = 1e-4)
  expect_equal(res$N, eaf_truth$N)


})

