library(dplyr)
test_that("multiplication works", {
  header_line <- "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ A1FREQ_CASES A1FREQ_CONTROLS INFO N N_CASES N_CONTROLS TEST BETA SE CHISQ LOG10P EXTRA"
  column_names <- unlist(strsplit(header_line, " "))


  tbl <- dplyr::tibble(cname = column_names[-9]) |> dplyr::mutate(val = 0) |>
    tidyr::pivot_wider(names_from = cname, values_from = val)

  guess_names(tbl)
})



test_that("download from gwas catalog work", {
  tbl <- tibble(
    CHROM       = character(),
    GENPOS      = integer(),
    ID          = character(),
    ALLELE1     = character(),
    ALLELE0     = character(),
    A1FREQ      = numeric(),
    BETA        = numeric(),
    N_CASES     = integer(),
    N_CONTROLS  = integer(),
    SE          = numeric(),
    N           = integer(),
    INFO        = numeric(),
    GC_ZSCORE = numeric(),
    NO_MATCH = numeric()
  )

  expect_no_error(guess_names(tbl))


})





test_that("REGENIE format is detected and renamed correctly", {
  re_df <- tibble(
    CHROM       = character(),
    GENPOS      = integer(),
    ID          = character(),
    ALLELE1     = character(),
    ALLELE0     = character(),
    A1FREQ      = numeric(),
    BETA        = numeric(),
    N_CASES     = integer(),
    N_CONTROLS  = integer(),
    INFO        = numeric(),
    SE          = numeric(),
    N           = integer()
  )

  out <- guess_names(re_df)

  expect_named(
    out,
    c("CHR","POS","RSID","EffectAllele","OtherAllele",
      "EAF","B","CaseN","ControlN","INFO","SE","N"),
    ignore.order = FALSE
  )
})

# -------------------------------------------------------------------------
# 2. Dictionary-only rescue (no template fits)
# -------------------------------------------------------------------------
test_that("dictionary rescues unmapped but critical columns", {
  dict_df <- tibble(
    SCAFFOLD  = character(),  # CHR synonym
    POSITION  = integer(),    # POS synonym
    SNP       = character(),  # RSID synonym
    A1        = character(),  # EffectAllele synonym
    A2        = character(),  # OtherAllele synonym
    BETA      = numeric(),
    STDERR    = numeric(),    # SE synonym
    PVAL      = numeric()
  )

  out <- guess_names(dict_df)

  expect_true(all(
    c("CHR","POS","RSID","EffectAllele","OtherAllele","B","SE","P") %in%
      names(out)
  ))
})


test_that("template plus dictionary yields complete set", {
  mixed_df <- tibble(
    CHROM      = character(),  # template
    GENPOS     = integer(),    # template
    ID         = character(),  # template
    ALLELE1    = character(),  # template
    ALLELE0    = character(),  # template
    EFFECT_BETA= numeric(),    # BETA synonym, *not* in template
    STDERR     = numeric(),    # SE synonym
    P_VALUE    = numeric(),    # P synonym
    rowid      = integer()     # should be dropped from matching
  )

  out <- guess_names(mixed_df)

  expect <- c("CHR","POS","RSID","EffectAllele","OtherAllele",
    "B","SE","P","rowid")

  expect_true(all(colnames(out) %in% expect))

})


test_that("unmapped columns persist with original names", {
  extra_df <- tibble(
    CHROM   = character(),
    GENPOS  = integer(),
    ID      = character(),
    ALLELE1 = character(),
    ALLELE0 = character(),
    BETA    = numeric(),
    SE      = numeric(),
    QCFLAG  = logical(),
    rowid   = integer()
  )

  out <- guess_names(extra_df)

  expect_true(all(c("QCFLAG","rowid") %in% names(out)))
  # and they are not renamed
  expect_identical(out$QCFLAG, logical())
  expect_identical(out$rowid, integer())
})
