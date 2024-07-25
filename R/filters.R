#' Remove rows with missing values, and write out the removed files to disk
#'
#' @param tbl a [dplyr::tibble()]
#' @param filepaths a list of filepaths, created by [setup_pipeline_paths()]
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' df <- remove_rows_with_na(sumstat, setup_pipeline_paths("testing"))
#' }
remove_rows_with_na <- function(tbl, filepaths) {


  # remove missing values from obligatory columns ---------------------------


  if(all(c("CHR", "POS") %in% colnames(tbl) & !"RSID" %in% colnames(tbl))) {

    tmp <- tidyr::drop_na(tbl, c("CHR", "POS", "EffectAllele", "OtherAllele"))

  } else if("RSID" %in% colnames(tbl) & !all(c("CHR", "POS") %in% colnames(tbl))) {

    tmp <- tidyr::drop_na(tbl, c("RSID", "EffectAllele", "OtherAllele"))

  } else if(all(c("CHR", "POS") %in% colnames(tbl) & "RSID" %in% colnames(tbl))) {
    # if we have all three, recode missing RSID to "." and it will be
    # updated if possible later
    tmp <- dplyr::mutate(tbl, RSID = dplyr::if_else(is.na(RSID), ".", RSID))

  }



  # identify removed rows ---------------------------------------------------
  na_rows <- dplyr::anti_join(tbl, tmp, by = "rowid") |>
    dplyr::select(rowid)

  if(nrow(na_rows) > 0) {

    outfile <- paste0(filepaths$removed, "missing_values.parquet")
    cli::cli_alert_danger("Found {nrow(na_rows)} rows with missing values. These are removed: ")
    cli::cli_inform("{.file {outfile}}")
    arrow::write_parquet(na_rows, outfile)

  } else {

    cli::cli_alert_success("No rows contained missing values")

  }

  tmp

}




#' Remove duplicated rows from a summary statistics file
#'
#' @param tbl a [dplyr::tibble()], with columns in `tidyGWAS()`
#' @param columns a character vector of columns passed to `dplyr::distinct()`
#' @param filepath a filepath to write out the removed rows
#'
#' @return a [dplyr::tibble()] with duplicates removed
#' @export
#'
#' @examples \dontrun{
#' remove_duplicates(tbl, c("CHR", "POS",), "duplicated_rows.tsv")
#' }

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
  cli::cli_alert_info("Lookins for duplications with columns: {columns}")
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





#' Detect Insertions/Deletions ('indels')
#' @description
#' Indels are detected by examining `EffectAllele` and `OtherAllele`
#'
#'
#' @param tbl a [dplyr::tibble()]
#' @inheritParams tidyGWAS
#' @inheritParams remove_rows_with_na
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' detect_indels(sumstat, TRUE, filepaths = setup_pipeline_paths("testing"))
#' }
detect_indels <- function(tbl, indel_strategy, filepaths, ...) {

  cli::cli_ol(c(
    "EffectAllele or OtherAllele, character length > 1: A vs AA",
    "EffectAllele or OtherAllele coded as 'D', 'I', or 'R'"
  ))


  tbl <- flag_indels(tbl)
  indels <-  dplyr::select(dplyr::filter(tbl,  .data[["indel"]]), -indel)
  tbl <- dplyr::select(dplyr::filter(tbl, !.data[["indel"]]), -indel)
  cli::cli_alert_success("Detected {nrow(indels)} rows as indels")

  if(nrow(indels) > 0) {

    indels <- validate_sumstat(
      tbl = indels,
      remove_cols = c("EffectAllele", "OtherAllele"),
      filter_func = make_callback(filepaths$removed_validate_indels),
      ...
    )
  }

  if(indel_strategy == "remove" & !is.null(indels)) {


    arrow::write_parquet(indels, filepaths$removed_indels)
    cli::cli_alert_warning("{.code keep_indels = FALSE}. Removed {nrow(indels)} rows as indels")
    cli::cli_inform("{.file {filepaths$removed_indels}}")
    indels <- NULL
  }

  list("main" = tbl, "indels" = indels)

}
