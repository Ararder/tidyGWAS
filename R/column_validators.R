utils::globalVariables(c("P", "p_was_0", "B", "EAF", "N", "missing_ea_oa", "SE", "non_acgt", "Z", "indel", "OR", "invalid_Z", "CaseN", "ControlN"))
impl_validators <- c("CHR", "POS", "EffectAllele", "OtherAllele","EAF", "SE", "P", "B", "Z", "N", "CaseN", "ControlN", "INFO")



# -------------------------------------------------------------------------

#' Validate statistics columns in a GWAS summary statistics file
#'
#' @param tbl a [dplyr::tibble()]
#' @param remove_cols Columns that should not be validated
#' @param filter_func handles reporting and writing removed files to disk
#' @inheritParams tidyGWAS
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' validate_sumstat(sumstat, remove_cols = "EffectAllele", convert_p = 0)
#' }
#'
validate_sumstat <- function(tbl, remove_cols = c(""), filter_func, convert_p) {


  stopifnot("remove_cols can only be a character vector" = is.character(remove_cols))

  if(!is.null(tbl)) {
    if(nrow(tbl) == 0) {
      tbl <- NULL
    }
  }


  # check that column validators and implemented, and which ones exist in tbl
  impl_validators <- impl_validators[!impl_validators %in% remove_cols]
  cols_in_sumstat <- colnames(tbl)[colnames(tbl) %in% impl_validators]

  for(colname in cols_in_sumstat) {
    tbl <- validate_columns(tbl = tbl, col = colname, convert_p = convert_p)
  }

  # use filter func if passed
  if(!missing(filter_func))  tbl <- filter_func(tbl)


  # finished ----------------------------------------------------------------

  tbl


}


#' Validate format of the RSID column in a GWAS summary statistics file
#'
#' @inheritParams tidyGWAS
#' @param filepath filepath to write out removed rows
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' validate_rsid(sumstat, "~/invalid_rsid.parquet")
#' }
validate_rsid <- function(tbl, filepath) {

  check_columns(c("RSID", "EffectAllele", "OtherAllele"), tbl)
  if(any(c("CHR", "POS") %in% colnames(tbl))) {
    stop("`validate_rsid should only be run on data.frames without CHR and POS")
  }

  tbl$RSID <- as.character(tbl$RSID)
  tbl <- flag_invalid_rsid(tbl)
  invalid_rsid_format <- dplyr::filter(tbl, invalid_rsid)


  # No invalid RSIDs --------------------------------------------------------

  if(nrow(invalid_rsid_format) == 0) {

    cli::cli_alert_success("All rows pass RSID validation")
    return(list("main" = tbl, "without_rsid" = NULL))

  }



  # -------------------------------------------------------------------------

  cli::cli_alert_info("Found {nrow(invalid_rsid_format)} rows with invalid RSID format: ")
  cli::cli_alert_info("Attempting to parse format...")


  # use EffectAllele and OtherAllele from sumstats, not from parsed format
  without_rsid <- split_rsid_by_regex(invalid_rsid_format) |>
    dplyr::select(-EffectAllele, -OtherAllele) |>
    dplyr::inner_join(dplyr::select(invalid_rsid_format, -RSID), by = "rowid")


  # print results to user, if any rows have been parsed correctly
  if(nrow(without_rsid) > 0) {
    cli::cli_alert_info("Parsed format")
    first_five_rows <- head(dplyr::select(without_rsid, RSID, CHR, POS, EffectAllele, OtherAllele), 5)
    cli::cat_print(first_five_rows, file = stderr())
  }


  # -------------------------------------------------------------------------


  failed <- dplyr::anti_join(invalid_rsid_format, without_rsid, by = "rowid")


  if(nrow(failed) > 0 ) {

    cli::cli_alert_danger("{nrow(failed)} rows had invalid RSID that could not be parsed and are removed.")
    cli::cli_inform("{.file {filepath}}")
    arrow::write_parquet(failed, filepath)

  }


  # clean up flags
  main <- dplyr::filter(tbl, !invalid_rsid) |> dplyr::select(-invalid_rsid)
  without_rsid <- dplyr::select(without_rsid, -RSID, -invalid_rsid)



  # failed parsed are removed

  list(
    "main" = main,
    "without_rsid" = without_rsid
  )

}





validate_columns <- function(tbl, col, convert_p) {

  if(col == "B") {

    tbl <- validate_b(tbl)

  } else if(col == "SE") {

    tbl <- validate_se(tbl)

  } else if(col == "EAF") {

    tbl <- validate_eaf(tbl)

  } else if(col == "N"){

    tbl <- validate_n(tbl)

  } else if(col == "CaseN"){

    tbl <- validate_casen(tbl)

  } else if(col == "ControlN"){

    tbl <- validate_controln(tbl)

  } else if(col == "Z") {

    tbl <- validate_z(tbl)

  } else if(col == "P") {

    tbl <- validate_p(tbl, convert_p = convert_p)


  } else if(col == "POS") {

    tbl <- validate_pos(tbl)

  } else if(col == "CHR") {
    tbl <- validate_chr(tbl)

  } else if(col == "EffectAllele" | col == "OtherAllele") {

    tbl <- validate_alleles(tbl, col = col)
  } else if(col == "INFO") {
    tbl <- validate_info(tbl)
  }

  # finished ----------------------------------------------------------------
  end_message(tbl, col = col)
  tbl

}




end_message <- function(tbl, col) {

  check <- !glue::glue("invalid_{col}") %in% colnames(tbl)
  if(check) {
    stop(cli::format_error(
      "Cannot print end msg for a column that has not been validated,
      invalid_{col} is not in {colnames(tbl)}
      ")
    )
  }

  n_invalid <- sum(tbl[[glue::glue("invalid_{col}")]])
  if(n_invalid > 0) {
    cli::cli_alert_warning("{n_invalid} rows failed {col} validation")
  } else {

  }

}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------



validate_chr <- function(tbl) {
  valid_chr <- c(1:22, "X", "Y", "MT", "XY")
  dplyr::mutate(
    tbl,
    CHR = as.character(CHR),
    CHR = stringr::str_to_upper(CHR),
    # can sometimes be chr22, or ch22
    CHR = stringr::str_remove(CHR, "CHR"),
    CHR = stringr::str_remove(CHR, "CH"),
    # can now handle
    CHR = dplyr::if_else(CHR == "XY", "XY", CHR),
    CHR = dplyr::if_else(CHR == "23", "X", CHR),
    CHR = dplyr::if_else(CHR == "M", "MT", CHR),
    invalid_CHR = dplyr::if_else(!CHR %in% valid_chr | is.na(CHR), TRUE, FALSE)
  )
}

validate_alleles <- function(tbl, col) {
  stopifnot(all(c("EffectAllele", "OtherAllele") %in% colnames(tbl)))
  stopifnot("EffectAllele and OtherAllele are not character columns, indicating mislabelling of columns" = is.character("EffectAllele") & is.character("OtherAllele"))
  possible_alleles <- c("A", "C","G","T")

  tbl <- dplyr::mutate(tbl, {{ col }} := stringr::str_to_upper(.data[[col]]))
  existing_alleles <- unique(tbl[[col]])
  if(!all(existing_alleles %in% possible_alleles)) cli::cli_alert_danger("Found non ACGT alleles in {col}: {existing_alleles}")
  new_col <- glue::glue("invalid_{col}")
  dplyr::mutate(tbl, {{ new_col }} := dplyr::if_else(!.data[[col]] %in% possible_alleles, TRUE, FALSE))
}

validate_pos <- function(tbl) {
  dplyr::mutate(
    tbl,
    POS = as.integer(POS),
    invalid_POS = dplyr::if_else(POS <= 0 | !is.finite(POS) | POS >= 10^9, TRUE, FALSE)
  )
}



validate_b <- function(tbl) {
  tbl$B <- as.double(tbl$B)
  median <- round(median(tbl$B, na.rm = TRUE), 5)
  if(dplyr::between(median, left = 0.9, right = 1.1)) cli::cli_alert_danger("WARNING: The median value of B is {median}, indicating that B has been mislabelled and contains odds-ratios")
  if(abs(median) > 0.1)  cli::cli_alert_danger("WARNING: The median value of B is {median}, which seems high")
  if(abs(median) <= 0.1)  cli::cli_alert_info("The median value of B is {median}, which seems reasonable")
  dplyr::mutate(tbl, invalid_B = dplyr::if_else(!is.finite(B), TRUE, FALSE))
}

validate_p <- function(tbl, convert_p) {
  dplyr::mutate(
    .data = tbl,
    P = as.double(P), p_was_0 = dplyr::if_else(P == 0, "Yes", "No"),
    P = dplyr::if_else(p_was_0 == "Yes", convert_p, P),
    invalid_P = dplyr::if_else(!is.finite(P) | P > 1 | P < 0, TRUE, FALSE)
  ) |>
    dplyr::select(-p_was_0)

}

validate_z <- function(tbl) {


  tbl <- dplyr::mutate(tbl, Z = as.double(Z))
  tbl <- dplyr::mutate(tbl, invalid_Z = dplyr::if_else(!is.finite(Z), TRUE, FALSE))
  maxval <- max(dplyr::filter(tbl, !invalid_Z)$Z)

  if(maxval < 100) {
    cli::cli_alert_info("Found {maxval} as largest absolute Z score, which seems reasonable")
  } else {
    cli::cli_alert_warning("WARNING: Found {maxval} as largest absolute Z score, which seems highy unlikely")
  }

  tbl
}


validate_se <- function(tbl) {
  dplyr::mutate(
    tbl,
    SE = as.double(SE),
    invalid_SE = dplyr::if_else(SE <= 0 | !is.finite(SE), TRUE, FALSE)
  )
}
validate_eaf <- function(tbl) {
  dplyr::mutate(
    tbl,
    EAF = as.double(EAF),
    invalid_EAF = dplyr::if_else(EAF <= 0 | EAF >= 1 | !is.finite(EAF), TRUE, FALSE)
  )
}

validate_n <- function(tbl) {
  dplyr::mutate(
    tbl,
    N = as.integer(N),
    invalid_N = dplyr::if_else(N <= 0 | !is.finite(N), TRUE, FALSE)
  )
}

validate_casen <- function(tbl) {
  dplyr::mutate(
    tbl,
    CaseN = as.integer(CaseN),
    invalid_CaseN = dplyr::if_else(CaseN <= 0 | !is.finite(CaseN), TRUE, FALSE)
  )

}

validate_controln <- function(tbl) {
  dplyr::mutate(
    tbl,
    ControlN = as.integer(ControlN),
    invalid_ControlN = dplyr::if_else(ControlN <= 0 | !is.finite(ControlN), TRUE, FALSE)
  )
}

validate_info <- function(tbl) {
  dplyr::mutate(
    tbl,
    INFO = as.double(INFO),
    # some leeway here, info is sometimes about 1
    invalid_INFO = dplyr::if_else(INFO <= 0 | INFO > 2 | !is.finite(INFO), TRUE, FALSE)
  )
}
