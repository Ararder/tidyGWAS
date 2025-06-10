test_that("repair_stats works can impute Z using B and P", {

  tmp <- dplyr::select(test_sumstat, -SE) |>
    dplyr::mutate(N = CaseN + ControlN)
  expect_no_error(repair_stats(tmp))

  # tmp2 <- dplyr::mutate(ll, SE = se_from_z_eaf_n(Z, EAF, N))
  # expect_equal(tmp2$SE, test_sumstat$SE, tolerance = 0.01)

})


test_that("If OR is provided and not B, OR should be converted to B", {
  tmp1 <- dplyr::mutate(pval_as_char_df, OR = exp(B)) |>
    dplyr::select(-B)
  repaired <- select_correct_columns(tmp1)


  expect_equal(repaired$B, pval_as_char_df$B)
})



test_that("If OR is provided and B exists, OR should be removed", {
  tmp1 <- dplyr::mutate(pval_as_char_df, OR = exp(B))
  repaired <- select_correct_columns(tmp1)
  expect_true(!"OR" %in% colnames(repaired))
})

test_that("Repair stats can repair B and SE from Z", {

  tmp <- dplyr::mutate(test_sumstat,N = CaseN + ControlN, Z = B/SE) |>
    dplyr::select(-B, -SE)
  after_repair <- repair_stats(tmp)

  cor(after_repair$B, test_sumstat$B)
  expect_true(cor(after_repair$B, test_sumstat$B) > 0.99)
  expect_true(cor(after_repair$SE, test_sumstat$SE) > 0.98)
  expect_true(cor(after_repair$Z, test_sumstat$B/test_sumstat$SE) > 0.99)


})


test_that("User should be informed if study N is used to impute B and SE", {


  z_only <- dplyr::mutate(test_sumstat, Z = B/SE) |>
    dplyr::select(-B, -SE)
  msg <- "*However, N does not not vary across SNPs, indicating it's the study-wide N and not SNP-wise N*"
  expect_message(regex = msg, repair_stats(dplyr::mutate(z_only,N  = 33433)))



})


test_that("repair_stats should be able to add P col if missing using Z", {
  test_sumstat <- flag_indels(pval_as_char_df) |> dplyr::filter(!indel)
  without_p <- dplyr::select(test_sumstat, -P)
  ll <- repair_stats(without_p)
  expect_equal(ll$P, as.double(test_sumstat$P), tolerance = 0.01)

})

test_that("repair_stats can add EAF", {


  test <- tidyGWAS::tidyGWAS(tbl, dbsnp_path) |>
    dplyr::select(CHR, POS_37, POS_38, RSID, EffectAllele, OtherAllele,SE)

  res <- repair_stats(tbl = test, dbsnp_path = dbsnp_path, impute_freq = "AFR")
  expect_true("EAF" %in% colnames(res))

  # can impute freq, and then calculate N
  expect_no_error(repair_stats(test, dbsnp_path, impute_freq = "EUR", impute_n = TRUE))

  # can impute N
  expect_no_error(repair_stats(tbl, impute_n = TRUE))


})



test_that("impute_n works", {
  file <- test_path("fixtures/HRC_eur_0.001.parquet")


  expect_no_error(res <- repair_stats(tbl, impute_n = TRUE))
  res <- repair_stats(dplyr::mutate(tbl, N = CaseN + ControlN), impute_n = TRUE)

  # mean(tbl$N), check that N is not different
  expect_equal(128002.1, mean(res$N), tolerance = 0.1)




})

