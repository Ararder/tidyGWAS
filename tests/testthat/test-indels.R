
test_that("detect indels work", {


  expect_no_error(res <- detect_indels(
    tbl = pval_as_char_df,
    indel_strategy = "qc",
    filepaths = filepaths,
    convert_p = convert_p,
    dbsnp_path = dbsnp_path
  ))

  expect_no_error(res <- detect_indels(
    tbl = tbl,
    indel_strategy = "qc",
    filepaths = filepaths,
    convert_p = convert_p,
    dbsnp_path = dbsnp_path
  ))



})




