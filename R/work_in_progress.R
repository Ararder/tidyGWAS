





validate_with_dbsnp <- function(data_list, build = c("NA", "37", "38"),filepaths, dbsnp_path) {

  rlang::check_required(dbsnp_path)
  build <- rlang::arg_match(build)
  cli::cli_h2("7a) Using dbSNP to repair and validate CHR:POS:RSID")
  cli::cli_li("Checking that EffectAllele and OtherAllele is compatible with REF and ALT in dbSNP")


  # -------------------------------------------------------------------------


  cli::cli_h3("Starting with main rows: ")
  data_list$main <- repair_dbnsp(data_list$main, dbsnp_path = dbsnp_path, build = build)
  filter_func <- make_callback(id = paste0(filepaths$removed_rows, "main_validate_with_dbsnp"))
  data_list$main <- filter_func(data_list$main)

  if(!is.null(data_list$without_rsid)) {
    cli::cli_h3("7b) rows without RSID: ")
    data_list$without_rsid <- repair_dbnsp(data_list$without_rsid, dbsnp_path = dbsnp_path, build = build)
    filter_func <- make_callback(id = "without_rsid_validate_with_dbsnp")
    data_list$without_rsid <- filter_func(data_list$without_rsid)

    # edge cases: -------------------------------------------------------------
    # it is possible that without_rsid subset contains rows that map to the same
    # rsid as already exists in main. Need to handle this

    data_list$main <- dplyr::bind_rows(data_list$main, data_list$without_rsid)
    data_list$without_rsid <- NULL
    before_unique_check <- dplyr::select(data_list$main, rowid)

    # handle duplications
    b38_missing <- dplyr::filter(data_list$main, is.na(CHR)) |>
      dplyr::distinct(CHR_37,POS_37,EffectAllele,OtherAllele, .keep_all = TRUE)

    data_list$main <- dplyr::filter(data_list$main, !is.na(CHR)) |>
      dplyr::distinct(CHR,POS,EffectAllele,OtherAllele, .keep_all = TRUE) |>
      dplyr::bind_rows(b38_missing)


    removed <- dplyr::anti_join(before_unique_check, data_list$main, by = "rowid")

    if(nrow(removed) > 0) {
      cli::cli_alert_info("Found {nrow(removed)} rows which map to a CHR:POS:REF:ALT that another variant maps to. These are removed")
      arrow::write_parquet(removed, paste0(filepaths$removed_rows, "map_to_dbsnp_duplications.parquet"))
    }


  }

  dplyr::bind_rows(data_list$main, data_list$indels)


}




repair_dbnsp <- function(tbl, dbsnp_path, build) {
  # existence of chr:pos or rsid decides which columns to repair
  has_rsid <- "RSID" %in% colnames(tbl)
  has_chr_pos <- all(c("CHR", "POS") %in% colnames(tbl))

  if(has_rsid & !has_chr_pos) {
    cli::cli_inform("Using RSID to align with dbSNP")
    tbl <- repair_chr_pos(tbl, dbsnp_path = dbsnp_path)

  } else if(has_chr_pos & !has_rsid) {
    cli::cli_inform("Using CHR and POS to align with dbSNP")
    tbl <- repair_rsid(tbl, build = build, dbsnp_path =  dbsnp_path)


  } else if(has_chr_pos & has_rsid) {
    cli::cli_inform("Using CHR, POS and RSID to align with dbSNP")
    tbl <- verify_chr_pos_rsid(tbl, build = build, dbsnp_path)

  }

  tbl
}
