

test_that("Repair_chr_pos works", {

  mock_dbsnp()
  test_file$rowid <- 1:nrow(test_file)
  tmp <- dplyr::select(test_file, -CHR, -POS)

  expect_no_error(repaired <- repair_chr_pos(tmp))

  # check that positions are imputed correctly
  in_dbsnp <- dplyr::filter(repaired, !no_dbsnp_entry) |>
    dplyr::select(CHR_37, POS_37, RSID, rowid)
  check <- dplyr::filter(test_file, rowid %in% in_dbsnp$rowid) |>
    dplyr::select(CHR, POS, RSID, rowid)
  merged <- dplyr::inner_join(in_dbsnp, check, by ="rowid") |>
    dplyr::mutate(CHR = as.character(CHR))
  expect_equal(merged$POS_37, merged$POS)
  expect_equal(merged$CHR_37, merged$CHR)

})


test_that("Repair_rsid works", {
  mock_dbsnp()
  tmp <-
    dplyr::mutate(test_file,  CHR = as.character(CHR), rowid = 1:nrow(test_file)) |>
    dplyr::select(rowid, CHR, POS, RSID, EffectAllele, OtherAllele)


  repaired <- repair_rsid(dplyr::select(tmp, -RSID))


  in_dbsnp <- dplyr::filter(repaired, !no_dbsnp_entry) |>
    dplyr::select(RSID, rowid)
  check <- dplyr::filter(tmp, rowid %in% in_dbsnp$rowid) |>
    dplyr::select(RSID, rowid)
  merged <- dplyr::inner_join(in_dbsnp, check, by ="rowid")
  expect_equal(merged$RSID.x, merged$RSID.y)




})



test_that("verify_chr_pos_rsid works", {
  mock_dbsnp()
  tmp <-
    dplyr::mutate(test_file,  CHR = as.character(CHR), rowid = 1:nrow(test_file)) |>
    dplyr::select(rowid, CHR, POS, RSID, EffectAllele, OtherAllele)

  expect_no_error(verify_chr_pos_rsid(tmp, bsgenome_objects = get_bsgenome()))


})





test_that("get_ref_data runs", {
  tmp_file <- withr::local_tempfile()
  data.table::fwrite(dplyr::tibble(RSID = "rs2020", new_RSID = "rs1090"), tmp_file)
  withr::local_envvar(c("rs_merge_arch" = tmp_file))
  expect_no_error(get_ref_data())

})

test_that("flatten_dbsnp runs", {
  expect_no_error(flatten_dbsnp(b38))
})

test_that("qc_with_dbsnp runs", {
  tmp <- dplyr::mutate(test_file, rowid = 1:nrow(test_file))
  expect_no_error(qc_with_dbsnp(tmp, b38))
})

test_that("make callback runs", {
  expect_no_error(make_callback(withr::local_tempfile()))
})
