test_that("multiplication works", {
  skip()

  false_id <- "GCST99999999"
  other_id <- "GCST90468178"
  third_id <- "GCST90245992"
  gwas_catalog_fpi("GCST90441306")
  gwas_catalog_fpi("GCST90468178")
  gwas_catalog_fpi(third_id)

})
