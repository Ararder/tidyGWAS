utils::globalVariables(c("discrep_freq", "base_EAF"))

remove_note <- function(x) {
  R.utils::dataFrame()
}

identify_removed_rows <- function(finished, filepaths) {
  # read in file before any munging was done
  start <- arrow::read_parquet(filepaths$raw_sumstats, col_select = "rowid")

  # find all removed rows
  removed_rows <- dplyr::anti_join(start, finished, by = "rowid")

  # get filepaths for all removed files
  files_in_dir <- list.files(
    paste0(filepaths$base, "/pipeline_info/"),
    pattern = "*removed_*",
    full.names = TRUE
  )

  removed_rows_with_flags <-
    files_in_dir |>
    purrr::set_names(base::basename) |>
    purrr::map(arrow::read_parquet) |>
    purrr::map(\(x) dplyr::select(x, rowid)) |>
    purrr::list_rbind(names_to = "reason") |>
    dplyr::select(rowid, reason) |>
    dplyr::mutate(reason = stringr::str_remove(reason, "removed_")) |>
    dplyr::mutate(reason = stringr::str_remove(reason, ".parquet"))

  if (sum(!removed_rows$rowid %in% removed_rows_with_flags$rowid)) {
    cli::cli_alert_danger(
      "WARNING: Could not track why some rows were removed. This means that some rows were removed
      that has not been accounted for. possibly unintended rows were removed."
    )
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

explain_removed_rows <- function(filepaths) {
  possible_reasons <- stringr::str_subset(filepaths, "removed") |>
    base::basename() |>
    stringr::str_remove("removed_") |>
    stringr::str_remove(".parquet")

  if (reason == "missing_alleles") {
    cli::cli_inform("Missing values in `EffectAllele` or `OtherAllele`")
  } else if (reason == "indels") {
    cli::cli_inform(
      "`indel_strategy='keep'` was not set resulting in indels being removed."
    )
  } else if (reason == "validate_indels") {
    cli::cli_inform(
      "Failed column validation, within the subset of rows detected as indels"
    )
  } else if (reason == "missing_critical") {
    cli::cli_inform("Missing values in either `CHR` or `POS`, or `RSID`")
  } else if (reason == "duplicates") {
    cli::cli_inform("Identified as duplicated rows")
  } else if (reason == "invalid_rsid") {
    cli::cli_inform(
      "Failed validation of RSID, and CHR and POS were not in the summary statistics"
    )
  } else if (reason == "validate_rsid_path") {
    cli::cli_inform("Failed column validation, only RSID detected")
  } else if (reason == "without_rsid") {}
}


write_finished_tidyGWAS <- function(df, output_format, filepaths) {
  df <- standardize_column_order(df)

  if (output_format == "hivestyle") {
    arrow::write_dataset(dplyr::group_by(df, CHR), filepaths$cleaned)
  } else if (output_format == "parquet") {
    arrow::write_parquet(
      df,
      paste(filepaths$base, "tidyGWAS_cleaned.parquet", sep = "/")
    )
  } else if (output_format == "csv") {
    arrow::write_csv_arrow(
      df,
      paste(filepaths$base, "tidyGWAS_cleaned.csv", sep = "/")
    )
  }
}

tidyGWAS_new_columns <- c(
  "POS_38",
  "POS_37",
  "REF_38",
  "REF_37",
  "indel",
  "multi_allelic",
  "rowid"
)
standardize_column_order <- function(tbl) {
  char_cols <- c("CHR", "EffectAllele", "OtherAllele")
  integer_cols <- c(
    "N",
    "CaseN",
    "ControlN",
    "EffectiveN",
    "POS_37",
    "POS_38",
    "rowid"
  )
  double <- c("B", "P", "EAF", "Z", "SE", "INFO")
  logical <- c("multi_allelic", "indel")

  first <- c("CHR", "POS_38", "POS_37", "RSID", "EffectAllele", "OtherAllele")
  second <- c(
    "B",
    "SE",
    "P",
    "EAF",
    "N",
    "CaseN",
    "ControlN",
    "EffectiveN",
    "Z",
    "INFO",
    "indel",
    "discrep_freq"
  )
  end <- c("rowid", "multi_allelic", "REF_37", "REF_38")

  dplyr::select(
    tbl,
    dplyr::all_of(first),
    dplyr::any_of(second),
    dplyr::all_of(end),
  ) |>
    dplyr::mutate(
      dplyr::across(dplyr::any_of(char_cols), ~ as.character(.)),
      dplyr::across(dplyr::any_of(integer_cols), ~ as.integer(.)),
      dplyr::across(dplyr::any_of(double), ~ as.double(.)),
      dplyr::across(dplyr::any_of(logical), ~ as.logical(.))
    )
}


check_columns <- function(columns, df) {
  missing <- columns[!columns %in% colnames(df)]
  name_missing <- names(columns[!columns %in% colnames(df)])
  if (!rlang::is_empty(missing)) {
    cli::cli_abort(
      "Provided column name are missing from data.frame:
      {.arg {missing}} = {.arg {name_missing}}.
      {.arg {missing}} is not in {colnames(df)}
      "
    )
  }
}


make_callback <- function(id) {
  outpath <- paste0(id, ".parquet")

  callback <- function(tbl) {
    # split into filter flags
    flags <- dplyr::select(tbl, rowid, dplyr::where(is.logical))
    # if ncol == 1, only rowid exists - no flags to filter on.
    if (ncol(flags) == 1) {
      cli::cli_inform("Found no flags to filter on")
      return(tbl)
    }

    remove <- dplyr::filter(
      flags,
      dplyr::if_any(dplyr::where(is.logical), \(x) x)
    )
    count_by_flag <-
      purrr::map(dplyr::select(remove, -rowid), \(x) sum(x, na.rm = T)) |>
      purrr::keep(\(x) x > 0)

    if (nrow(remove) > 0) {
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


update_rsid <- function(tbl, filepath, dbsnp_path) {
  # detect merged rsIDs -----------------------------------------------------

  dset <- arrow::open_dataset(paste(
    dbsnp_path,
    "refsnp-merged/part-0.parquet",
    sep = "/"
  ))
  updates <- dplyr::select(tbl, dplyr::all_of(c("rowid", "RSID"))) |>
    dplyr::semi_join(dset, y = _, by = c("old_RSID" = "RSID")) |>
    dplyr::collect()

  rsid_info <-
    dplyr::left_join(tbl, updates, by = c("RSID" = "old_RSID")) |>
    dplyr::mutate(
      new_RSID = dplyr::if_else(!is.na(RSID.y), RSID.y, RSID),
      old_RSID = dplyr::if_else(!is.na(RSID.y), RSID, NA_character_)
    ) |>
    dplyr::select(rowid, RSID = new_RSID, old_RSID, -RSID.y)

  # add updated rsid to tbl
  tbl <- dplyr::inner_join(
    dplyr::select(tbl, -RSID),
    rsid_info,
    by = "rowid"
  ) |>
    dplyr::select(-old_RSID)

  # identify rows with updated rsid
  updated_rows <- sum(!is.na(rsid_info$old_RSID))

  if (updated_rows > 0) {
    cli::cli_alert_success("{updated_rows} rows with updated RSID")
    cli::cli_li("{.file {filepath}}")
    arrow::write_parquet(dplyr::filter(rsid_info, !is.na(old_RSID)), filepath)
  } else {
    cli::cli_li("Found no RSIDs that has been merged")
  }

  tbl
}


select_correct_columns <- function(tbl) {
  # check input columns
  cli::cli_h3("Checking that columns follow tidyGWAS format")
  cli::cli_alert_success(
    "The following columns are used for further steps: {.emph {colnames(tbl)[colnames(tbl) %in% valid_column_names]}}"
  )

  # check for invalid columns
  n_invalid_cols <- length(
    colnames(tbl)[!colnames(tbl) %in% valid_column_names] > 0
  )
  if (n_invalid_cols) {
    cli::cli_alert_danger(
      "{.strong Removed columns:  {colnames(tbl)[!colnames(tbl) %in% valid_column_names]}}"
    )
  }
  tbl <- dplyr::select(tbl, dplyr::any_of(valid_column_names))

  # remove columns with all NA
  cli::cli_h3("Checking for columns with all NA")
  na_cols <- colnames(tbl)[colSums(is.na(tbl)) == nrow(tbl)]
  tbl <- dplyr::select(tbl, -dplyr::all_of(na_cols))

  if (length(na_cols) > 0) {
    cli::cli_alert_danger(
      "The following columns were removed as they contained only NA's:
                         {.emph {na_cols}}"
    )
  } else {
    cli::cli_alert_success("Found no columns with all NA")
  }

  # To continue, we need at least one of the following sets of columns
  no_chr_pos <- !all(c("CHR", "POS") %in% colnames(tbl))
  missing_rsid <- !"RSID" %in% colnames(tbl)
  if (no_chr_pos & missing_rsid) {
    stop("Either CHR and POS or RSID are required columns")
  }
  if (!all(c("EffectAllele", "OtherAllele") %in% colnames(tbl))) {
    stop("EffectAllele and OtherAllele are required columns")
  }

  # handle B/OR -------------------------------------------------------------

  if ("OR" %in% colnames(tbl) & !"B" %in% colnames(tbl)) {
    cli::cli_alert_info(
      "Found OR but not BETA. converting to B using {.code base::log(OR)}"
    )
    tbl <- dplyr::mutate(tbl, B = log(OR)) |>
      dplyr::select(-OR)
  }

  if ("OR" %in% colnames(tbl) & "B" %in% colnames(tbl)) {
    cli::cli_alert_info("found OR and B, removing OR")
    tbl <- dplyr::select(tbl, -OR)
  }

  # handle N ----------------------------------------------------------------

  if (
    all(c("CaseN", "ControlN") %in% colnames(tbl)) &
      !"EffectiveN" %in% colnames(tbl)
  ) {
    cli::cli_alert_info(
      "Found CaseN and ControlN, and no effective N:
                        Calculating EffectiveN by {.code EffectiveN = 4 / (1 / ControlN + 1 / CaseN)}"
    )
    round(tbl$EffectiveN <- 4 / (1 / tbl$ControlN + 1 / tbl$CaseN))
  }

  if (
    all(c("CaseN", "ControlN") %in% colnames(tbl)) & !"N" %in% colnames(tbl)
  ) {
    cli::cli_alert_info(
      "Found CaseN and ControlN, Calculating N by {.code N = ControlN + CaseN}"
    )
    tbl$N <- tbl$CaseN + tbl$ControlN
  }

  if (!"N" %in% colnames(tbl)) {
    cli::cli_alert_danger(
      "Found no N column. It is highly recommended to supply a value for N, as many downstream GWAS applications rely on this information"
    )
  }

  tbl
}



apply_dbsnp_filter <- function(tbl, filepaths) {
  n_before <- nrow(tbl)
  before_filters <- tbl
  tbl <- dplyr::filter(tbl, !no_dbsnp_entry & !incompat_alleles) |>
    dplyr::select(-dplyr::all_of(c("no_dbsnp_entry", "incompat_alleles")))

  removed_no_dbsnp <- dplyr::anti_join(before_filters, tbl, by = "rowid")
  if (nrow(removed_no_dbsnp) > 0) {
    cli::cli_alert_warning(
      "Removed {nrow(removed_no_dbsnp)} rows with no dbSNP entry or with incompat alleles"
    )
    cli::cli_inform("{.file {filepaths$removed_no_dbsnp}}")
    arrow::write_parquet(removed_no_dbsnp, filepaths$removed_no_dbsnp)
  }

  tbl
}


read_freq_ref <- function(flag_discrep_freq, dbsnp_path) {

  df_eaf <- arrow::open_dataset(paste0(dbsnp_path, "/EAF_REF_1KG")) |>
    dplyr::filter(.data[["ancestry"]] == flag_discrep_freq) |>
    dplyr::select(CHR, POS_38 = POS, EffectAllele, OtherAllele, EAF) |>
    dplyr::collect()



}

add_freq_diff_flag <- function(main, flag_discrep_freq, dbsnp_path) {

  df_eaf <- read_freq_ref(flag_discrep_freq,dbsnp_path) |>
    dplyr::distinct()


  sub <- dplyr::select(main, CHR, POS_38, EffectAllele, OtherAllele, base_EAF = EAF)





  full_sub <-
    dplyr::inner_join(sub, df_eaf, by = c("CHR","POS_38", "EffectAllele" ="OtherAllele", "OtherAllele" = "EffectAllele")) |>
    dplyr::mutate(EAF = 1-EAF) |>
    dplyr::select(c("CHR","POS_38", "EffectAllele", "OtherAllele", "base_EAF", "EAF")) |>
    dplyr::bind_rows(dplyr::inner_join(sub, df_eaf, by = c("CHR","POS_38", "EffectAllele", "OtherAllele"))) |>
    dplyr::mutate(
      diff = abs(base_EAF - EAF),
      discrep_freq = dplyr::if_else(diff >= 0.2, TRUE, FALSE)
    ) |>
    dplyr::select(CHR,POS_38, EffectAllele, OtherAllele, discrep_freq)

  final <- dplyr::left_join(main, full_sub, by = c("CHR","POS_38", "EffectAllele", "OtherAllele"))

  counts <- dplyr::count(final, .data[["discrep_freq"]])

  cli::cli_alert_success(
    "Matched {sum(counts[1,2]$n, counts[2,2]$n)}/{nrow(final)} variants with reference data allele frequency."
    )

  cli::cli_alert_info(
    "{round(counts[2,2]$n / nrow(final),1)*100}% ({counts[2,2]$n}) variants with absolute allele frequency difference of >= 20%"
  )

  final
}

#' #' Save the filepath for the dbSNP reference data
#' #'
#' #' @param dbsnp_path filepath to dbSNP
#' #'
#' #' @returns NULL
#' #' @export
#' #'
#' #' @examples \dontrun{
#' #' set_default_dbsnp_path("data/local/dbSNP155")
#' #' }
#' set_default_dbsnp_path <- function(dbsnp_path) {
#'   rlang::is_scalar_character(dbsnp_path) || cli::cli_abort("dbsnp_path must be a character string.")
#'   file.exists(dbsnp_path) || cli::cli_abort("dbsnp_path must be a valid file path.")
#'
#'   target <- file.path(Sys.getenv("HOME"), ".config", "dbSNP155")
#'
#'   if(!file.exists(dirname(target))) dir.create(dirname(target))
#'
#'   file.symlink(from = path.expand(dbsnp_path), to = target)
#'   cli::cli_alert_success(
#'     "Set dbSNP155 path to {.file {dbsnp_path}}.
#'     You can change this by running {.code set_default_dbsnp_path()} again."
#'   )
#' }
