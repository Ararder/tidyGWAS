remove_rows_with_na <- function(tbl, columns, filepath) {

  tmp <- tidyr::drop_na(tbl, dplyr::all_of(columns))
  na_rows <- dplyr::anti_join(tbl, tmp, by = "rowid") |>
    dplyr::select(rowid)

  if(nrow(na_rows) > 0) {


    cli::cli_alert_danger("Found {nrow(na_rows)} rows with missing values in {columns}. These are removed: ")
    cli::cli_inform("{.file {filepath}}")
    arrow::write_parquet(na_rows, filepath)

  } else {

    cli::cli_alert_success("No rows contained missing values in {columns}")

  }

  tmp

}






remove_duplicates <- function(tbl, columns = NULL, filepath) {


  # arrange by P if possible - dplyr::distinct() keeps the first row
  if("P" %in% colnames(tbl)) tbl <- dplyr::arrange(tbl, .data[["P"]])

  # -------------------------------------------------------------------------

  if(is.null(columns)) {
    if(all(c("CHR", "POS") %in% colnames(tbl))) {
      columns <- c("CHR", "POS", "EffectAllele", "OtherAllele")
    } else {
      columns <- c("RSID", "EffectAllele", "OtherAllele")
    }
  }

  # -------------------------------------------------------------------------



  # remove duplicates using cols_to_use
  no_dups <- dplyr::distinct(tbl, dplyr::pick(dplyr::all_of(columns)), .keep_all = TRUE)

  # find which rows were removed
  removed <- dplyr::anti_join(tbl, no_dups, by = "rowid")

  # logging
  cli::cli_alert_info("Looking for duplications with columns: {columns}")
  if(nrow(removed) > 0) {

    cli::cli_alert_danger("Removed {nrow(removed)} rows flagged as duplications")
    cli::cli_li("{.file {filepath}}")
    arrow::write_parquet(removed, filepath)

  } else {

    cli::cli_alert_success("Found no duplications")
  }

  # done --------------------------------------------------------------------


  no_dups

}





detect_indels <- function(tbl, indel_strategy, filepaths, convert_p) {

  cli::cli_ol(c(
    "EffectAllele or OtherAllele, character length > 1: A vs AA",
    "EffectAllele or OtherAllele coded as 'D', 'I', or 'R'"
  ))


  tbl <- flag_indels(tbl)
  indels <-  dplyr::select(dplyr::filter(tbl,  .data[["indel"]]), -indel)
  without_indels <- dplyr::select(dplyr::filter(tbl, !.data[["indel"]]), -indel)
  cli::cli_alert_success("Detected {nrow(indels)} rows as indels")


  if(indel_strategy == "remove") {

    if(!is.null(indels)) {
      arrow::write_parquet(indels, filepaths$removed_indels)
      cli::cli_alert_warning("{.code indel_strategy = 'remove'}. Removed {nrow(indels)} rows as indels")
      cli::cli_inform("{.file {filepaths$removed_indels}}")
    }
    # yes this is intentional, to return an empty tibble of same format
    indels <- dplyr::filter(indels, FALSE)
  } else if(nrow(indels) > 0 & indel_strategy == "keep") {

    indels <- validate_sumstat(
      tbl = indels,
      remove_cols = c("EffectAllele", "OtherAllele"),
      filter_func = make_callback(filepaths$removed_validate_indels),
      convert_p = convert_p
    )

  }


  list("main" = without_indels, "indels" = indels)

}
