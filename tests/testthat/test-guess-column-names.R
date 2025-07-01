test_that("multiplication works", {
  header_line <- "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ A1FREQ_CASES A1FREQ_CONTROLS INFO N N_CASES N_CONTROLS TEST BETA SE CHISQ LOG10P EXTRA"
  column_names <- unlist(strsplit(header_line, " "))


  tbl <- dplyr::tibble(cname = column_names[-9]) |> dplyr::mutate(val = 0) |>
    tidyr::pivot_wider(names_from = cname, values_from = val)

  guess_names(tbl)
})



test_that("download from gwas catalog work", {
  skip()
  dbnsp <- "~/Downloads/dbSNP155/"
  path <- "/var/folders/6f/ndpsm1z96js7qcvmf6vtz638stnpv6/T/RtmpAmgHXl/GCST90475332/GCST90475332.tsv.gz"

  test <- tidyGWAS(
    tbl = path,
    delim = "\t",
    dbsnp = dbnsp
  )

})
