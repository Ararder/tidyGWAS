#' Remove all columns that do not follow tidyGWAS naming
#'
#' @param tbl a [dplyr::tibble()] as created by [parse_tbl()]
#' @return a [dplyr::tibble()], with only columns following [tidyGWAS_columns()] naming kept
#' @export
#'
#' @examples \dontrun{
#' sumstats <- select_correct_columns(sumstats)
#' }
select_correct_columns <- function(tbl) {

  # check input columns
  cli::cli_h3("Checking that columns follow tidyGWAS format")
  cli::cli_alert_success("The following columns are used for further steps: {.emph {colnames(tbl)[colnames(tbl) %in% valid_column_names]}}")

  n_invalid_cols <- length(colnames(tbl)[!colnames(tbl) %in% valid_column_names] > 0)
  if(n_invalid_cols) cli::cli_alert_danger("{.strong Removed columns:  {colnames(tbl)[!colnames(tbl) %in% valid_column_names]}}")

  tbl <- dplyr::select(tbl, dplyr::any_of(valid_column_names))

  # In some formats, columns are kept with all NA, to provide a consistent structure
  # of the file. For tidyGWAS, this is not necessary, and we remove these columns.
  cli::cli_h3("Checking for columns with all NA")
  na_cols <- colnames(tbl)[colSums(is.na(tbl)) == nrow(tbl)]
  tbl <- dplyr::select(tbl, -dplyr::all_of(na_cols))

  if(length(na_cols) > 0 ){
  cli::cli_alert_danger("The following columns were removed as they contained only NA's:
                         {.emph {na_cols}}")
  } else {
    cli::cli_alert_success("Found no columns with all NA")
  }

  # To continue, we need at least one of the following sets of columns
  if(all(!c("CHR", "POS") %in% colnames(tbl)) & !"RSID" %in% colnames(tbl)) stop("Either CHR and POS or RSID are required columns")
  if(!all(c("EffectAllele", "OtherAllele") %in% colnames(tbl))) stop("EffectAllele and OtherAllele are required columns")



  # handle B/OR -------------------------------------------------------------

  if("OR" %in% colnames(tbl) & !"B" %in% colnames(tbl)) {
    cli::cli_alert_info("Found OR but not BETA. converting to B using {.code base::log(OR)}")
    tbl <- dplyr::mutate(tbl, B = log(OR)) |>
      dplyr::select(-OR)
  }

  if("OR" %in% colnames(tbl) & "B" %in% colnames(tbl)) {
    cli::cli_alert_info("found OR and B, removing OR")
    tbl <- dplyr::select(tbl, -OR)
  }


  # handle N ----------------------------------------------------------------

  if(all(c("CaseN", "ControlN") %in% colnames(tbl))) tbl$N <- (tbl$CaseN + tbl$ControlN)


  if(!"N" %in% colnames(tbl)) cli::cli_alert_danger("Found no N column, and no study_n was supplied. It is highly recommended to supply a value for N, as many downstream GWAS applications rely on this information")


  tbl

}


#' Remove rows with missing values, and write out the removed files to disk
#'
#' @param tbl a [dplyr::tibble()] as created by [parse_tbl()]
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

#' Remove duplicated rows
#' @description
#' remove_duplicates uses either CHR:POS:EffectAllele:OtherAllele
#' or RSID:EffectAllele:OtherAllele to compute uniqueness.
#'
#' If possible rows are arranged by p-value, to select the row with the smallest P.
#'
#'
#' @inheritParams remove_rows_with_na
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' paths <- setup_pipeline_paths("testing")
#' df <- remove_duplicates(sumstat, paths)
#' }
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


    cli::cli_alert_danger("Removed {nrow(removed)} rows flagged as duplications")
    cli::cli_li("{.file {filepaths$removed_duplicates}}")
    arrow::write_parquet(removed, filepaths$removed_duplicates)

  } else {

    cli::cli_alert_success("Found no duplications")

  }

  no_dups
}

#' Update rsIDs from dbSNP that have been merged into other RSIDs
#'
#' @inheritParams remove_rows_with_na
#' @param dbsnp_path filepath to dbSNP155 directory
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' update_rsid(sumstat, filepaths = setup_pipeline_paths("testing"), dbsnp_path = "~/dbSNP155")
#' }
update_rsid <- function(tbl, filepaths, dbsnp_path) {


  # detect merged rsIDs -----------------------------------------------------

  dset <- arrow::open_dataset(paste(dbsnp_path, "refsnp-merged", sep = "/"))
  updates <- dplyr::select(tbl, dplyr::all_of(c("rowid", "RSID"))) |>
    dplyr::semi_join(dset, y = _,  by = c("old_RSID" = "RSID")) |>
    dplyr::collect()

  rsid_info <-
    dplyr::left_join(tbl, updates, by = c("RSID" = "old_RSID")) |>
    dplyr::mutate(
      new_RSID = dplyr::if_else(!is.na(RSID.y), RSID.y, RSID),
      old_RSID = dplyr::if_else(!is.na(RSID.y), RSID, NA_character_)
    ) |>
    dplyr::select(rowid, RSID = new_RSID, old_RSID, -RSID.y)

  # add updated rsid to tbl
  tbl <- dplyr::inner_join(dplyr::select(tbl, -RSID), rsid_info, by = "rowid") |>
    dplyr::select(-old_RSID)

  # identify rows with updated rsid
  updated_rows <- sum(!is.na(rsid_info$old_RSID))

  if(updated_rows > 0) {

    cli::cli_alert_success("{updated_rows} rows with updated RSID")
    cli::cli_li("{.file {filepaths$updated_rsid}}")
    arrow::write_parquet(dplyr::filter(rsid_info, !is.na(old_RSID)), filepaths$updated_rsid)

  } else {

    cli::cli_li("Found no RSIDs that has been merged")

  }

  tbl

}




#' Detect Insertions/Deletions ('indels')
#' @description
#' Indels are detected by examining `EffectAllele` and `OtherAllele`
#'
#'
#' @param tbl a [dplyr::tibble()] as created by [parse_tbl()]
#' @inheritParams tidyGWAS
#' @inheritParams remove_rows_with_na
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' detect_indels(sumstat, TRUE, filepaths = setup_pipeline_paths("testing"))
#' }
detect_indels <- function(tbl, indel_strategy, filepaths,...) {

  cli::cli_ol(c(
    "EffectAllele or OtherAllele, character length > 1: A vs AA",
    "EffectAllele or OtherAllele coded as 'D', 'I', or 'R'"
  ))


  tbl <- flag_indels(tbl)
  indels <-  dplyr::select(dplyr::filter(tbl,  .data[["indel"]]), -indel)
  tbl <- dplyr::select(dplyr::filter(tbl, !.data[["indel"]]), -indel)
  cli::cli_alert_success("Detected {nrow(indels)} rows as indels")

  indels <- validate_sumstat(
    tbl = indels,
    remove_cols = c("EffectAllele", "OtherAllele"),
    filter_func = make_callback(filepaths$removed_validate_indels),
    id = "indel_rows",
    ...
  )

  if(indel_strategy == "remove" & !is.null(indels)) {


    arrow::write_parquet(indels, filepaths$removed_indels)
    cli::cli_alert_warning("{.code keep_indels = FALSE}. Removed {nrow(indels)} rows as indels")
    cli::cli_inform("{.file {filepaths$removed_indels}}")
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
      cli::cli_h3("Listing how many rows are removed per flag: ")
      cli::cli_dl(purrr::list_simplify(count_by_flag))
      cli::cli_li("Removed a total of {nrow(remove)} rows: {.file {outpath}}")
    } else {
      cli::cli_h3("{.emph All rows passed validation}")
    }



    arrow::write_parquet(remove, outpath)

    dplyr::select(tbl, -dplyr::where(is.logical)) |>
      dplyr::filter(!rowid %in% remove$rowid)
  }

  callback
}

apply_filters <- function(tbl, filepaths) {


  # no dbSNP mapping
  clean <- dplyr::filter(tbl, !no_dbsnp_entry)

  removed_no_dbsnp <- dplyr::anti_join(tbl, clean, by = "rowid")

  # missing on either build
  clean2 <- tidyr::drop_na(clean, c("CHR", "POS", "CHR_37", "POS_37", "RSID"))
  removed_missing_on_either_build <- dplyr::anti_join(clean, clean2, by = "rowid")

  # chr mismatch between builds
  clean <- dplyr::filter(clean2, CHR_37 == CHR)
  removed_chr_mismatch <- dplyr::anti_join(clean2, clean, by = "rowid")


  # communicate removed rows ------------------------------------------------

  if(nrow(removed_no_dbsnp) > 0) {
    cli::cli_alert_warning("Removed {nrow(removed_no_dbsnp)} rows with no dbSNP entry")
    cli::cli_inform("{.file {filepaths$removed_no_dbsnp}}")
    arrow::write_parquet(removed_no_dbsnp, filepaths$removed_no_dbsnp)
  }

  if(nrow(removed_missing_on_either_build) > 0) {
    cli::cli_alert_warning("Removed {nrow(removed_missing_on_either_build)} rows with missing on either build")
    cli::cli_inform("{.file {filepaths$removed_missing_on_either_build}}")
    arrow::write_parquet(removed_missing_on_either_build, filepaths$removed_missing_on_either_build)
  }

  if(nrow(removed_chr_mismatch) > 0) {
    cli::cli_alert_warning("Removed {nrow(removed_chr_mismatch)} rows with chr mismatch between builds")
    cli::cli_inform("{.file {filepaths$removed_rows_chr_mismatch}}")
    arrow::write_parquet(removed_chr_mismatch, filepaths$removed_chr_mismatch)
  }


  # -------------------------------------------------------------------------



  clean |>
    dplyr::select(-incompat_alleles, -no_dbsnp_entry)

}
