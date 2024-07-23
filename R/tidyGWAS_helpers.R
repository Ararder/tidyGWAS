parse_tbl <- function(tbl, ...) {


  if("character" %in% class(tbl)) {
    stopifnot("A filepath has to be a character vector of length 1" = "character" %in% class(tbl) & length(tbl) == 1)
    stopifnot("File does not exist"  = file.exists(tbl))
    filename <- basename(tbl)
    md5 <- tools::md5sum(tbl)
    tbl <- arrow::read_delim_arrow(tbl, ...)

    ncols <- length(colnames(tbl))
    if(ncols <= 3) {
      stop(glue::glue("Only {ncols} columns in the file. At least 3 columns are required.
                      Most likely you forgot to use the `delim` argument to specify the delimiter."))
    }


  } else if("data.frame" %in% class(tbl)) {
    md5 <- NULL
    tbl <- dplyr::tibble(tbl)
    filename <- "raw"

  } else {

    stop("tbl is not a character vector or a data.frame")
  }

  if(!"rowid" %in% colnames(tbl)) tbl$rowid <- 1:nrow(tbl)

  list("tbl" = tbl, "filename" = filename, "md5" = md5)
}


identify_removed_rows <- function(finished, filepaths) {

  start <- arrow::read_parquet(filepaths$raw_sumstats, col_select = "rowid")

  removed_rows <- dplyr::anti_join(start, finished, by = "rowid")

  # remove the file with all rowids
  files_in_dir <- list.files(paste0(filepaths$base, "/pipeline_info/"), pattern = "*removed_*", full.names = TRUE)


  removed_rows_with_flags <-
    files_in_dir |>
    purrr::set_names(base::basename) |>
    purrr::map(arrow::read_parquet) |>
    purrr::map(\(x) dplyr::select(x, rowid)) |>
    purrr::list_rbind(names_to = "reason") |>
    dplyr::select(rowid, reason) |>
    dplyr::mutate(reason = stringr::str_remove(reason, "removed_")) |>
    dplyr::mutate(reason = stringr::str_remove(reason, ".parquet"))

  if(sum(!removed_rows$rowid %in% removed_rows_with_flags$rowid)) {
    cli::cli_alert_danger(
      "WARNING: Could not track why some rows were removed. This means that some rows were removed
      that has not been accounted for. possibly unintended rows were removed.")
  }

  breakdown <-
    removed_rows_with_flags |>
    dplyr::semi_join(removed_rows, by = "rowid") |>
    dplyr::count(reason)

  vec <- c(breakdown$n)
  names(vec) <- breakdown$reason
  cli::cli_h3("Listing final breakdown of removed rows: ")
  cli::cli_dl(vec)




}

write_finished_tidyGWAS <- function(df, output_format, filepaths) {

  df <- standardize_column_order(df)

  if(output_format == "hivestyle") {

    arrow::write_dataset(dplyr::group_by(df, CHR), filepaths$cleaned)

  } else if(output_format == "parquet") {

    arrow::write_parquet(df, paste(filepaths$base, "tidyGWAS_cleaned.parquet", sep = "/"))

  } else if(output_format == "csv") {

    arrow::write_csv_arrow(df, paste(filepaths$base, "tidyGWAS_cleaned.csv", sep = "/"))

  }


}


standardize_column_order <- function(tbl) {
  dplyr::select(tbl, dplyr::any_of(
    c(
      "CHR", "POS", "RSID", "EffectAllele", "OtherAllele", "EAF",
      "Z", "B", "SE", "P", "N", "CaseN", "ControlN", "INFO"
      ,"CHR_37", "POS_37", "rowid", "multi_allelic", "indel",
      "REF" = "ref_allele")), dplyr::everything()
  )

}

update_column_names <- function(tbl, column_map, CaseN = NULL, ControlN =NULL, N=NULL) {
  rlang::check_required(tbl)


  if(!missing(column_map)) {
    stopifnot(
      "All entries in column_map must be present in the input tbl" =
        all(column_map %in% colnames(tbl))
    )
    tbl <- dplyr::rename(tbl, !!!column_map)
  }

  if(!is.null(CaseN)) tbl <- dplyr::mutate(tbl, CaseN =  {{ CaseN }})
  if(!is.null(ControlN)) tbl <- dplyr::mutate(tbl, ControlN =  {{ ControlN }})
  if(!is.null(N)) tbl <- dplyr::mutate(tbl, N = {{ N }})



  tbl
}


create_id <- function(tbl, build = c("37", "38")) {
  build <- rlang::arg_match(build)
  ref_name <- paste0("REF_", build)
  tbl |>
    dplyr::mutate(
      ref_allele = .data[[ref_name]],
      POS = .data[[paste0("POS_", build)]]
    ) |> 
    dplyr::mutate(
      alt_allele = dplyr::if_else(ref_allele == EffectAllele, OtherAllele, EffectAllele),
      ID = stringr::str_c(CHR, POS, ref_allele, alt_allele, sep = ":")
    ) |>
    dplyr::select(-"alt_allele", -"ref_allele", -"POS")

}

