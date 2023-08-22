test_that("flag_indel works", {

  tmp <- flag_indels(pval_as_char_df) |>
    dplyr::filter(indel) |>
    nrow()
  expect_true(tmp == 1839)


})




test_that("flag_incorrect_rsid_format detects incorrect formats", {
  valid_rsids <-   c("rs2321", "rs12331230","Rs100", "rS94", "RS12313")
  not_valid_rsids <- c("ss12313", "X:123213:A:C", "Y_12123123_A_T", "14:043824", "rs12331230_A_C", "14:043824:rs12331230")
  examples <- rep(
    c(valid_rsids, not_valid_rsids),
    length.out = 27715
    )

  detected <- dplyr::mutate(pval_as_char_df, RSID = examples) |>
    flag_invalid_rsid() |>
    dplyr::filter(!invalid_rsid) |>
    dplyr::count(RSID) |>
    dplyr::pull(RSID)

  expect_true(
    all(detected %in% valid_rsids)
  )


})

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
  repaired <- repair_stats(tmp1)


  expect_equal(repaired$B, pval_as_char_df$B)
})



test_that("If OR is provided and B exists, OR should be removed", {
  tmp1 <- dplyr::mutate(pval_as_char_df, OR = exp(B))
  repaired <- repair_stats(tmp1)
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


test_that("split_rsid_by_regex correctly splits rows", {
  valid <-   c("Y:10050_T_A", "X:10050_T_A", "14_96458_C_G", "14:96458:C:G", "13-1400000")
  not_valid <- c("ss12313", "X:123XY13:A:C", "LL_12123123_A_T", "14:043824YXZ", "rs12331230_A_C", "14:043824:rs12331230")
  examples <- rep(
    c(valid, not_valid),
    length.out = 27715
  )

  matches <- dplyr::mutate(pval_as_char_df, RSID = examples, rowid = 1:nrow(pval_as_char_df)) |>
    split_rsid_by_regex() |>
    dplyr::count(RSID) |>
    dplyr::pull(RSID)

  # i expect all elements in valid to be in matches,
  # and none of the elements in not_valid to be in matches
  expect_true(all(valid %in% matches) & all(!not_valid %in% matches))

})



test_that("split_rsid_by_regex correctly splits rows", {

  res <- pval_as_char_df |>
    dplyr::mutate(
      expEA = dplyr::case_when(
      EffectAllele == "A" ~ "T",
      EffectAllele == "T" ~ "A",
      EffectAllele == "G" ~ "C",
      EffectAllele == "C" ~ "G",
      .default = EffectAllele
    ),
    expOA = dplyr::case_when(
      OtherAllele == "A" ~ "T",
      OtherAllele == "T" ~ "A",
      OtherAllele == "G" ~ "C",
      OtherAllele == "C" ~ "G",
      .default = OtherAllele
      )
    ) |>
    strand_flip() |>
    dplyr::mutate(
      t1 = EffectAllele == expEA,
      t2 = OtherAllele == expOA,
    )

  expect_true(all(res$t1 & res$t2))
})




test_that("flag duplicates work", {
  n_dup_rsid <- dplyr::add_row(pval_as_char_df, dplyr::slice(pval_as_char_df, 100)) |>
    flag_duplicates() |>
    dplyr::summarise(dup = sum(dup_rsid))

  # n_dup_chr_pos <-
  tmp <- dplyr::add_row(pval_as_char_df, dplyr::slice(pval_as_char_df, 100))

  tmp[27716, ]$RSID <- "rs7812154071"
  n_dup_chr_pos <-
    flag_duplicates(tmp, column = "chr_pos") |>
    dplyr::summarise(dup = sum(dup_chr_pos))

  expect_true(n_dup_chr_pos$dup == 2)
  expect_true(n_dup_rsid$dup == 2)

})

