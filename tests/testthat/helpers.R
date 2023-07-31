

randomizer <- function(end=10) {
  sample(c(1:10),nrow(tbl), replace=TRUE) == 1
}
mock_dbsnp <- function(){
  tbl <- tmp |>
    dplyr::slice_sample(n = 50000)


  save(tbl, file = test_path("fixtures/mri_test_data.rda"))
}


