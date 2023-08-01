
b38 <- map_to_dbsnp(dplyr::tibble(test_file), by = "rsid", build = 38)
b37 <- map_to_dbsnp(dplyr::tibble(test_file), by = "rsid", build = 37)


save(b38, file = test_path("fixtures/b38.rds"))
save(b37, file = test_path("fixtures/b37.rds"))