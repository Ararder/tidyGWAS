load(test_path("data/sumstats/test_sumstat.rds"))
load(test_path("data/sumstats/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))


# repair_stats ------------------------------------------------------------


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
  z_only <- dplyr::mutate(test_file,N = CaseN + ControlN, Z = B/SE) |>
    dplyr::select(-B, -SE)

  repaired <- repair_stats(z_only)
  expect_true("B" %in% colnames(repaired))
  expect_true("SE" %in% colnames(repaired))
  expect_equal(test_file$B, repaired$B, tolerance = 0.05)
  expect_equal(test_file$SE, repaired$SE, tolerance = 0.05)

})


test_that("User should be informed if study N is used to impute B and SE", {


  z_only <- dplyr::mutate(test_file, Z = B/SE) |>
    dplyr::select(-B, -SE)
  msg <- "*However, N does not not vary across SNPs, indicating it's the study-wide N and not SNP-wise N*"
  expect_message(regex = msg, repair_stats(dplyr::mutate(z_only,N  = 33433)))



})


test_that("repair_stats should be able to add P col if missing using Z", {
  without_p <- dplyr::select(test_file, -P)
  ll <- repair_stats(without_p)
  expect_equal(ll$P, test_file$P, tolerance = 0.005)

})


#  repair_chr_pos and repair_rsid ----------------------------------------



test_that("Repair_chr_pos works", {
  skip("Time consuming")
  bsgenome_objects <- get_bsgenome()
  testrep <- repair_chr_pos(dplyr::select(tbl, -CHR, -POS))



})


test_that("Repair_rsid works", {
  skip("Time consuming")
  bsgenome_objects <- get_bsgenome()
  tbl <- dplyr::mutate(tbl, CHR = as.character(CHR), POS = as.integer(POS))
  testrep <- repair_rsid(dplyr::select(tbl, -RSID), bsgenome_objects)

})



