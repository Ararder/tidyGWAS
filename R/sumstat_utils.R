utils::globalVariables(c("new_RSID", "old_RSID"))

#' Find all rows which are part of a set of duplicated rows
#' @description
#' Many duplication tools such as [base::duplicated()] or [dplyr::distinct()]
#' identify rows which are duplications. It is often useful to see ALL rows
#' which are part of the duplication set, and not just the second row.
#'
#' creates new column: `dup_rsid` or `dup_chr_pos`, a T/F flag.
#' Specifically, flags both rows in a duplication pair, and not just first or
#' last duplicate row, making it easy to work with all rows that are part of a
#' duplication
#'
#' @param tbl a [dplyr::tibble()]
#' @param column Which columns should be used to form a unique ID?
#'
#' @return a tibble with a new column marking duplicates
#' @export
#'
#' @examples \dontrun{
#'
#' # will tag multi-allelics as duplications
#' flag_duplicates(tbl, column = "rsid")
#' flag_duplicates(tbl, column = "chr_pos")
#' # if you are interested in rows that are variant duplications
#' flag_duplicates(tbl, column = "rsid_ref_alt")
#' flag_duplicates(tbl, column = "chr_pos_ref_alt")
#'
#' }
flag_duplicates <- function(tbl, column = c("rsid", "chr_pos", "chr_pos_ref_alt", "rsid_ref_alt")) {
  column = rlang::arg_match(column)
  new_name <- glue::glue("dup_{column}")

  if(column == "rsid") {

    columns <- c("RSID")
    dplyr::mutate(tbl, {{ new_name }} := (duplicated(tbl[,columns]) | duplicated(tbl[,columns], fromLast = TRUE)))

  } else if(column == "chr_pos") {

    columns <- c("CHR", "POS")
    dplyr::mutate(tbl, {{ new_name }} := (duplicated(tbl[,columns]) | duplicated(tbl[,columns], fromLast = TRUE)))

  } else if(column == "chr_pos_ref_alt") {

    columns <- c("CHR", "POS", "EffectAllele", "OtherAllele")
    dplyr::mutate(tbl, {{ new_name }} := (duplicated(tbl[,columns]) | duplicated(tbl[,columns], fromLast = TRUE)))

  } else if(column == 'rsid_ref_alt') {

    columns <- c("RSID", "EffectAllele", "OtherAllele")
    dplyr::mutate(tbl, {{ new_name }} := (duplicated(tbl[,columns]) | duplicated(tbl[,columns], fromLast = TRUE)))


  }
}

#' Detect "indels" in GWAS summary statistics
#'
#' @param tbl a [dplyr::tibble()] with columns `EffectAllele` and `OtherAllele`
#'
#' @return a [dplyr::tibble()] with a TRUE/FALSE column `indel` added, where
#' indel == TRUE corresponds to a row marked as an indel.
#' @export
#'
#' @examples \dontrun{
#' all_indels <-
#'   flag_indels(tbl) |>
#'   dplyr::filter(indels)
#' }
flag_indels <- function(tbl) {
  stopifnot(all(c("EffectAllele", "OtherAllele") %in% colnames(tbl)))
  indel_allele_codes <- c("D", "I", "R")

  tbl <- dplyr::mutate(tbl,
    EffectAllele = stringr::str_to_upper(EffectAllele),
    OtherAllele = stringr::str_to_upper(OtherAllele),
    indel = (stringr::str_length(.data[["EffectAllele"]]) + stringr::str_length(.data[["OtherAllele"]])) > 2,
    indel = dplyr::if_else(indel, indel, EffectAllele %in% indel_allele_codes | OtherAllele %in% indel_allele_codes)
    )

  tbl

}

#' Detect entries that are not valid rsID's in GWAS summary statistics
#'
#' @param tbl a [dplyr::tibble()] with column `RSID`.
#' @param regex regex used to detect non-RSIDs
#'
#' @return a [dplyr::tibble()] with column `invalid_rsid`
#' @export
#'
#' @examples \dontrun{
#' flag_invalid_rsid(tbl) |>
#' dplyr::filter(invalid_rsid)
#' }
flag_invalid_rsid <- function(tbl, regex = "^[rR][sS]?\\d{1,10}$") {
  stopifnot("RSID" %in% colnames(tbl))
  dplyr::mutate(tbl, invalid_rsid = stringr::str_detect(.data[["RSID"]], regex, negate=TRUE))

}


split_rsid_by_regex <- function(tbl) {
  stopifnot(all(c("RSID", "rowid") %in% colnames(tbl)))

  df <- dplyr::select(tbl, dplyr::all_of(c("rowid", "RSID")))

  sep_regex <- "[:_-]"
  regex <- list(
    "^(?:[1-9]|1[0-9]|2[0-2]|X|Y)[:_-][1-9][0-9]{0,9}[:_-][A-Za-z]{1,100}[:_-][A-Za-z]{1,100}$",
    "^(?:[1-9]|1[0-9]|2[0-2]|X|Y)[:_-][1-9][0-9]{0,9}$",
    "^chr(?:[1-9]|1[0-9]|2[0-2]|[XxYy])[:_-][1-9][0-9]{0,9}$"
  )

  colnames <- list(
    c("CHR","POS", "EffectAllele", "OtherAllele"),
    c("CHR","POS"),
    c("CHR","POS")
  )



  out <- purrr::map2(regex, colnames, \(reg, column_names)
              # only attempt to separate wider when rows match regex
              dplyr::filter(df ,stringr::str_detect(.data[["RSID"]], reg)) |>
                tidyr::separate_wider_delim(
                  cols = RSID,
                  delim = stringr::regex(sep_regex),
                  names = column_names,
                  cols_remove = FALSE
                )
              ) |>
  # remove empty lists and rowbind results
  purrr::keep(\(x) nrow(x) > 0) |>
  purrr::list_rbind() |>
  dplyr::tibble()

  cli::cli_alert_success("Succeded in parsing {nrow(out)} out of {nrow(tbl)} rows")

  # Need to provide a consistent out format.
  if(!"EffectAllele" %in% colnames(out)) out$EffectAllele <- NA_character_
  if(!"OtherAllele" %in% colnames(out)) out$OtherAllele <- NA_character_
  if(!"CHR" %in% colnames(out)) out$CHR <- NA_character_
  if(!"POS" %in% colnames(out)) out$POS <- NA_integer_
  if(!"rowid" %in% colnames(out)) out$rowid <- NA_integer_
  if(!"RSID" %in% colnames(out)) out$RSID <- NA_character_

  out

}




#' Strand flip alleles
#'
#' @param tbl a [dplyr::tibble()] with columns `EffectAllele` and `OtherAllele`
#'
#' @return a [dplyr::tibble()] with columns `EffectAllele` and `OtherAllele` flipped
#' @export
#'
#' @examples \dontrun{
#' tbl <- strand_flip(tbl)
#' }
strand_flip <- function(tbl) {
  stopifnot("Requires columns EffectAllele & OtherAllele" =
              all(c("EffectAllele", "OtherAllele") %in% colnames(tbl))
  )
  dplyr::mutate(tbl,
                EffectAllele = dplyr::case_when(
                  EffectAllele == "A" ~ "T",
                  EffectAllele == "T" ~ "A",
                  EffectAllele == "G" ~ "C",
                  EffectAllele == "C" ~ "G",
                  .default = EffectAllele
                ),
                OtherAllele = dplyr::case_when(
                  OtherAllele == "A" ~ "T",
                  OtherAllele == "T" ~ "A",
                  OtherAllele == "G" ~ "C",
                  OtherAllele == "C" ~ "G",
                  .default = OtherAllele
                )
  )
}



