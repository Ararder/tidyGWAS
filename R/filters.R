select_correct_columns <- function(tbl, study_n) {
  # check input columns
  cli::cli_h3("1) Checking that columns follow tidyGWAS format")
  cli::cli_alert_success("The following columns are used for further steps: {.emph {colnames(tbl)[colnames(tbl) %in% valid_column_names]}}")
  n_invalid_cols <- length(colnames(tbl)[!colnames(tbl) %in% valid_column_names] > 0)
  if(n_invalid_cols) cli::cli_alert_danger("{.strong Removed columns:  {colnames(tbl)[!colnames(tbl) %in% valid_column_names]}}")

  tbl <- dplyr::mutate(tbl, rowid = as.integer(row.names(tbl)))
  tbl <- dplyr::select(tbl,dplyr::any_of(valid_column_names))

  if(all(!c("CHR", "POS") %in% colnames(tbl)) & !"RSID" %in% colnames(tbl)) {
    stop("Either CHR and POS or RSID are required columns")
  }
  if(!all(c("EffectAllele", "OtherAllele") %in% colnames(tbl))) {
    stop("EffectAllele and OtherAllele are required columns")
  }

  # handle N ----------------------------------------------------------------

  if(all(c("CaseN", "ControlN") %in% colnames(tbl))) tbl$N <- (tbl$CaseN + tbl$ControlN)

  if(!missing(study_n))  {
    cli::cli_alert_info("Using N = {study_n} as N")
    tbl <- dplyr::mutate(tbl, N = {{ study_n }})
  }

  if(!"N" %in% colnames(tbl)) cli::cli_alert_danger("Found no N column, and no study_n was supplied. It is highly recommended to supply a value for N, as many downstream GWAS applications rely on this information")


  tbl

}


remove_rows_with_na <- function(tbl, filepaths) {


  tmp <- tidyr::drop_na(tbl)
  na_rows <- dplyr::anti_join(tbl, tmp, by = "rowid") |>
    dplyr::select(rowid)

  if(nrow(na_rows) > 0) {
    outfile <- paste0(filepaths$removed, "missing_values.parquet")
    cli::cli_alert_danger("Found {nrow(na_rows)} rows with missing values. These are removed: ")
    cli::cli_inform("{.file {outfile}}")
    arrow::write_parquet(na_rows, outfile)
  }

  tmp

}

remove_duplicates <- function(tbl, filepaths) {


  # arrange by P if possible - so that row with smallest pval is select
  # in case of duplicate
  if("P" %in% colnames(tbl)) tbl <- dplyr::arrange(tbl, .data[["P"]])

  if(all(c("CHR", "POS") %in% colnames(tbl))) {
    id <- "CHR_POS_REF_ALT"
    cols_to_use <- c("CHR", "POS", "EffectAllele", "OtherAllele")
  } else {
    id <- "RSID_REF_ALT"
    cols_to_use <- c("RSID", "EffectAllele", "OtherAllele")
  }

  # remove duplicates using cols_to_use
  no_dups <- dplyr::distinct(tbl, dplyr::pick(dplyr::all_of(cols_to_use)), .keep_all = TRUE)

  # find which rows were removed
  removed <- dplyr::anti_join(tbl, no_dups, by = "rowid")



  cli::cli_alert_info("A unique ID is formed by concontenating {cols_to_use}")
  if(nrow(removed) > 0) {
    out <- paste(filepaths$removed_rows, "duplicated_rows.parquet")
    cli::cli_alert_danger("Removed {nrow(removed)} rows flagged as duplications")
    cli::cli_li("{.file {out}}")
    arrow::write_parquet(removed, out, compression = "gzip")
  }

  no_dups
}

update_rsid <- function(tbl, filepaths, dbsnp_path) {

  rsid_info <- update_merged_rsid(tbl, paste(dbsnp_path, "utils","RsMergeArch.parquet", sep = "/"))
  tbl <- dplyr::inner_join(dplyr::select(tbl, -RSID), rsid_info, by = "rowid") |>
    dplyr::select(-old_RSID)

  cli::cli_li("{sum(!is.na(rsid_info$old_RSID))} rows with updated RSID: {.file {filepaths$updated_rsid}}")
  arrow::write_parquet(dplyr::filter(rsid_info, !is.na(old_RSID)), filepaths$updated_rsid)

  tbl
}



detect_indels <- function(tbl, keep_indels, filepaths) {

  cli::cli_ol(c(
    "EffectAllele or OtherAllele, character length > 1: A vs AA",
    "EffectAllele or OtherAllele coded as 'D', 'I', or 'R'"
  ))


  tbl <- flag_indels(tbl)
  indels <-  dplyr::select(dplyr::filter(tbl,  .data[["indel"]]), -indel)
  tbl <- dplyr::select(dplyr::filter(tbl, !.data[["indel"]]), -indel)
  if(nrow(indels) == 0) indels <- NULL
  if(!isTRUE(keep_indels) & !is.null(indels)) {

    outpath <- paste0(filepaths$removed_rows, "indels_removed.parquet")
    arrow::write_parquet(indels, outpath)
    cli::cli_alert_warning("{.code keep_indels = FALSE}. Removed {nrow(indels)} rows as (indels)")
    cli::cli_inform("{.file {outpath}}")
    indels <- NULL
  }

  list("main" = tbl, "indels" = indels)

}


make_callback <- function(id) {


  outpath <- paste0(id, ".parquet")

  callback <- function(tbl) {
    # split into filter flags
    flags <- dplyr::select(tbl, rowid, dplyr::where(is.logical))
    # if ncol == 1, only rowid exists - no flags to filter on.
    if(ncol(flags) == 1) {
      cli::cli_inform("Found no flags to filter on")
      return(tbl)
    }

    remove <- dplyr::filter(flags, dplyr::if_any(dplyr::where(is.logical), \(x) x))
    count_by_flag <-
      purrr::map(dplyr::select(remove, -rowid), \(x) sum(x, na.rm = T)) |>
      purrr::keep(\(x) x > 0)



    if(nrow(remove) > 0) {
      cli::cli_h2("Listing how many rows are removed per flag: ")
      cli::cli_dl(purrr::list_simplify(count_by_flag))
      cli::cli_inform("Removed a total of {nrow(remove)} rows: {.file {outpath}}")
    } else {
      cli::cli_h2("{.emph All rows passed validation}")
    }



    arrow::write_parquet(remove, outpath, compression = "gzip")

    dplyr::select(tbl, -dplyr::where(is.logical)) |>
      dplyr::filter(!rowid %in% remove$rowid)
  }

  callback
}
