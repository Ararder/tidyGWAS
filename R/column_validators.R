utils::globalVariables(c("P", "p_was_0", "B", "EAF", "N", "missing_ea_oa", "SE", "non_acgt", "Z", "indel", "OR"))
impl_validators <- c("CHR", "POS", "EffectAllele", "OtherAllele","EAF", "SE", "P", "B", "Z", "N")



# -------------------------------------------------------------------------

#' Validate statistics columns in a GWAS summary statistics file
#'
#' @param tbl a [dplyr::tibble()] as created by [parse_tbl()]
#' @param remove_cols Columns that should not be validated
#' @param filter_func handles reporting and writing removed files to disk
#' @inheritParams tidyGWAS
#' @param id Used to customize messages.
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' validate_sumstat(sumstat, remove_cols = "EffectAllele", convert_p = 0)
#' }
#'
validate_sumstat <- function(tbl, remove_cols= c(""), filter_func,  verbose = FALSE, convert_p, id) {

  # check if 0 rows tbl or NULL
  stopifnot("remove_cols can only be a character vector" = is.character(remove_cols))
  tbl <- check_zero_rows(tbl)
  if(is.null(tbl)) return(NULL)


  # check that column validators and implemented, and which ones exist in tbl
  impl_validators <- impl_validators[!impl_validators %in% remove_cols]
  cols_in_sumstat <- colnames(tbl)[colnames(tbl) %in% impl_validators]

  for(colname in cols_in_sumstat) {
    tbl <- validate_columns(tbl = tbl, col = colname, verbose = verbose, convert_p = convert_p)
  }

  # use filter func if passed
  if(!missing(filter_func))  tbl <- filter_func(tbl)


  # finished ----------------------------------------------------------------

  tbl


}


#' Validate format of the RSID column in a GWAS summary statistics file
#'
#' @param tbl a [dplyr::tibble()] as created by [parse_tbl()]
#' @inheritParams tidyGWAS
#' @param outpath Filepath: Where to write rows with invalid RSID?
#'
#' @return a tbl
#' @export
#'
#' @examples \dontrun{
#' validate_rsid(sumstat, "~/invalid_rsid.parquet")
#' }
validate_rsid <- function(tbl, verbose = FALSE, outpath) {

  if(verbose) start_message("RSID")
  tbl$RSID <- as.character(tbl$RSID)
  tbl <- flag_invalid_rsid(tbl)
  invalid_rsid_format <- dplyr::filter(tbl, invalid_rsid)


  # No invalid RSIDs --------------------------------------------------------

  if(nrow(invalid_rsid_format) == 0) {
    cli::cli_alert_success("All rows pass RSID validation")
    return(list("main" = tbl, "without_rsid" = NULL))
  }


  # invalid RSIDs, but other columns exist ----------------------------------

  if(all(c("CHR", "POS", "EffectAllele", "OtherAllele") %in% colnames(tbl))) {
    cli::cli_alert_info("{ nrow(dplyr::filter(tbl, invalid_rsid)) } rows had an invalid RSID. RSID will be repaired using dbSNP if possible")
    return(
      list(
        "main" = dplyr::filter(tbl, !invalid_rsid),
        "without_rsid" = dplyr::filter(tbl, invalid_rsid) |> dplyr::select(-RSID, -invalid_rsid)
      ))
  }


  # need to parse RSID ------------------------------------------------------
  # no CHR or POS exists in tbl

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
    cli::cli_inform("{.file {outpath}}")
    arrow::write_parquet(failed, outpath)

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







#' Check that values in GWAS summary statistics columns are correct
#'
#' `validate_columns()` does not remove any rows, but adds a TRUE/FALSE
#' flag for the specified column.
#'
#'
#'
#'
#' @param tbl a [dplyr::tibble()] as created by [parse_tbl()]
#' @param col Which column to check values in?
#' @inheritParams tidyGWAS
#'
#'
#' @return a [dplyr::tibble()], with a column added named as invalid_{col}.
#' If you validate the "B" column, validate_columns will add a TRUE/FALSE column
#' named invalid_B to the input tibble.
#' @export
#' @examples \dontrun{
#' gwas_file <- validate_columns(
#'  tbl = gwas_file,
#'  col = "B",
#'  verbose = FALSE,
#'  # if you want to keep 0 pvalues as 0.
#'  convert_p = 0
#' )
#' dplyr::filter(gwas_file, invalid_P)
#' }
#'
validate_columns <- function(
    tbl,
    col = c("B", "SE", "EAF", "N", "Z", "P","POS","CHR", "EffectAllele", "OtherAllele"),
    verbose = TRUE,
    convert_p = 2.225074e-308
    ) {
  col = rlang::arg_match(col)
  if(verbose) start_message(col)

  if(col == "B") {

    tbl$B <- as.double(tbl$B)
    median <- round(median(tbl$B, na.rm = TRUE), 5)
    if(dplyr::between(median, left = 0.9, right = 1.1)) cli::cli_alert_danger("WARNING: The median value of B is {median}, indicating that B has been mislabelled and contains odds-ratios")
    if(abs(median) > 0.1)  cli::cli_alert_danger("WARNING: The median value of B is {median}, which seems high")
    if(abs(median) <= 0.1)  cli::cli_alert_info("The median value of B is {median}, which seems reasonable")

    tbl <- dplyr::mutate(tbl, invalid_B = dplyr::if_else(!is.finite(B), TRUE, FALSE))

  } else if(col == "SE") {

    tbl <- dplyr::mutate(tbl, SE = as.double(SE), invalid_SE = dplyr::if_else(SE <= 0 | !is.finite(SE), TRUE, FALSE))

  } else if(col == "EAF") {

    tbl <-dplyr::mutate(tbl, EAF = as.double(EAF), invalid_EAF = dplyr::if_else(EAF <= 0 | EAF >= 1 | !is.finite(EAF), TRUE, FALSE))

  } else if(col == "N"){

    tbl <- dplyr::mutate(tbl,N = as.integer(N),invalid_N = dplyr::if_else(N <= 0 | !is.finite(N), TRUE, FALSE))

  } else if(col == "Z") {

    tbl <- dplyr::mutate(tbl, Z = as.double(Z))
    maxval <- max(abs(tbl$Z))

    if(maxval < 100) {
      cli::cli_alert_info("Found {maxval} as largest absolute Z score, which seems reasonable")
    } else {
      cli::cli_alert_warning("WARNING: Found {maxval} as largest absolute Z score, which seems highy unlikely")
    }
    tbl <- dplyr::mutate(tbl, invalid_Z = dplyr::if_else(!is.finite(Z), TRUE, FALSE))

  } else if(col == "P") {

    tbl <- dplyr::mutate(
      .data = tbl,
      P = as.double(P), p_was_0 = dplyr::if_else(P == 0, "Yes", "No"),
      P = dplyr::if_else(p_was_0 == "Yes", convert_p, P),
      invalid_P = dplyr::if_else(!is.finite(P) | P > 1 | P < 0, TRUE, FALSE)
      ) |>
      dplyr::select(-p_was_0)



  } else if(col == "POS") {

    tbl <- dplyr::mutate(tbl, POS = as.integer(POS), invalid_POS = dplyr::if_else(POS <= 0 | !is.finite(POS) | POS >= 10^9, TRUE, FALSE))

  } else if(col == "CHR") {

    valid_chr <- c(1:22, "X", "Y", "MT", "XY")
    tbl <- dplyr::mutate(tbl,
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

  } else if(col == "EffectAllele" | col == "OtherAllele") {

    stopifnot(all(c("EffectAllele", "OtherAllele") %in% colnames(tbl)))
    stopifnot("EffectAllele and OtherAllele are not character columns, indicating mislabelling of columns" = is.character("EffectAllele") & is.character("OtherAllele"))
    possible_alleles <- c("A", "C","G","T")

    tbl <- dplyr::mutate(tbl, {{ col }} := stringr::str_to_upper(.data[[col]]))
    existing_alleles <- unique(tbl[[col]])
    if(!all(existing_alleles %in% possible_alleles)) cli::cli_alert_danger("Found non ACGT alleles in {col}: {existing_alleles}")
    new_col <- glue::glue("invalid_{col}")
    tbl <- dplyr::mutate(tbl, {{ new_col }} := dplyr::if_else(!.data[[col]] %in% possible_alleles, TRUE, FALSE))
    tbl
  }

  # finished ----------------------------------------------------------------
  end_message(tbl, col = col)
  tbl

}






start_message <- function(col) {
  if(col == "P") {

    cli::cli_h3("Validating the P column:")
    cli::cli_ol()
    cli::cli_li("Coerce to double {.code base::as.double()}")
    cli::cli_li("If P == 0, P is converted to {.arg convert_p}")
    cli::cli_li("checks for NA, NaN or Inf")
    cli::cli_li("checks that P is within the range: P >= 0 & P > 1")
  }

 if(col == "B") {
   cli::cli_h3("Validating the B column:")
   cli::cli_ol()
   cli::cli_li("Will coerce to double: base::as.double()")
   cli::cli_li("Check for mislabelled OR by looking at median value")
   cli::cli_li("checks for NA, NaN or Inf")
   cli::cli_li("Print median value")
 }

  if(col == "SE") {
    cli::cli_h3("Validating the SE column:")
    cli::cli_ol()
    cli::cli_li("Will coerce to double: base::as.double()")
    cli::cli_li("Check for SE <= 0")
    cli::cli_li("Checks for NA, NaN or Inf")
  }

  if(col == "EAF") {
    cli::cli_h3("Validating the EAF column:")
    cli::cli_ol()
    cli::cli_li("Will coerce to double: base::as.double()")
    cli::cli_li("Check for EAF <= 0 & EAF >= 1")
    cli::cli_li("Checks for NA, NaN or Inf")
  }
  if(col == "N") {
    cli::cli_h3("Validating the N column:")
    cli::cli_ol()
    cli::cli_li("Will coerce to integer: base::as.integer()")
    cli::cli_li("Check for N <= 0")
    cli::cli_li("Checks for NA, NaN or Inf")
  }
  if(col == "Z") {
    cli::cli_h3("Validating the Z column:")
    cli::cli_ol()
    cli::cli_li("Will coerce to double: base::as.integer()")
    cli::cli_li("Check for max(abs(z) >= 100")
    cli::cli_li("Checks for NA, NaN or Inf")
  }
  if(col == "RSID") {
    cli::cli_h3("Validating the RSID column:")
    cli::cli_ul()
    cli::cli_li("Checking that RSID follows rs format: [rR][sS][1-10]")
    cli::cli_li("If rows fail rs format, look for CHR:POS or CHR:POS:REF:ALT format")
  }

  if(col == "EffectAllele" | col == "OtherAllele") {
    cli::cli_h3("Validating {col}")
    cli::cli_ol()
    cli::cli_li("Will error if not type = character")
    cli::cli_li("{col} converted to uppercase")
    cli::cli_li("Check for values that are not A,C,G or T")
    cli::cli_li("Check for NA")

  }

  if(col == "CHR") {
    cli::cli_h3("Validating CHR")
    cli::cli_ol()
    cli::cli_li("CHR coerced to character {.code base::as.character()}")
    cli::cli_li("converted to uppercase {.code stringr::str_to_upper()}")
    cli::cli_li("'CHR' removed (handles UCSC style CHR format)")
    cli::cli_li("Check for values outside 1:22, X, Y")
    cli::cli_li("Check for NA")

  }

  if(col == "POS") {
    cli::cli_h3("Validating POS")
    cli::cli_ol()
    cli::cli_li("POS coerced to integer {.code base::as.integer()}")
    cli::cli_li("Checks for NA, NaN, Inf")
    cli::cli_li("Checks for POS <= 0")
    cli::cli_li("Check for POS >= 10^9")


  }



}

end_message <- function(tbl, col) {

  stopifnot("Cannot print end msg for column that doesnt exist" = glue::glue("invalid_{col}") %in% colnames(tbl))
  n_invalid <- sum(tbl[[glue::glue("invalid_{col}")]])
  if(n_invalid > 0) {
    cli::cli_alert_warning("{n_invalid} rows failed {col} validation")
  } else {
    # cli::cli_alert_success("All rows pass {col} validation")

  }

}
