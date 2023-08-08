utils::globalVariables(c(
  "new_RSID", "old_RSID", "retracted","reactivated", ":="
))

flag_indels <- function(tbl) {
  stopifnot(all(c("EffectAllele", "OtherAllele") %in% colnames(tbl)))
  indel_allele_codes <- c("D", "I", "R")

  tbl <- dplyr::mutate(tbl,
    EffectAllele = stringr::str_to_upper(EffectAllele),
    OtherAllele = stringr::str_to_upper(OtherAllele),
    indel = (stringr::str_length(.data[["EffectAllele"]]) + stringr::str_length(.data[["OtherAllele"]])) > 2,
    indel = dplyr::if_else(indel, indel, EffectAllele %in% indel_allele_codes | OtherAllele %in% indel_allele_codes)
    )

  cli::cli_alert_success("Detected {sum(tbl$indel)} rows as indels")
  tbl

}
flag_incorrect_rsid_format <- function(tbl, regex = "^[rR][sS]?\\d{1,10}$") {
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
# some reference on what these files are:
# SNPhistory columns: https://www.biostars.org/p/157447/
# What is going on with these files? https://www.ncbi.nlm.nih.gov/books/NBK573473/

# downloaded from NCBI 2023-07-04, dbsnp151 (why not dbsnp155? i dont understand why they are not available there...)
# without_rs_grch38:
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606/known_issues/rs_without_GRCh38_mapping_b151.bcp.README
# SNPhistory
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/database/organism_data/SNPHistory.bcp.gz
# RsMerge Arch
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/database/organism_data/RsMergeArch.bcp.gz

# TEMP DURING DEVELOPMENT: THESE FILES ARE HERE ON LONGLEAF
# rs_merge_arch_path = "~/arvhar/extdata/RsMergeArch.parquet",
# snp_history_path="~/arvhar/extdata/SNPhistory.parquet.gz",
# without_rs_grch38="~/arvhar/extdata/b151_rs_without_GRCh38_mapping.parquet.gz",

#' Repair statistics column in a GWAS summary statistics tibble
#'
#' `tidyGWAS::repar_stats()` is a collection of functions that can be used to
#' infer missing columns in GWAS summary statistics. The functions are based on
#' functionality found online.
#'
#' `repair_stats()` uses the tidyGWAS column specification. The functionality
#' has been tested assuming the input data.frame is a tibble.
#' @param tbl input a tibble (data frame)
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' updated <- repair_stats(my_gwas)
#' }
repair_stats <- function(tbl) {
  # reparation of Z from B and P is from LDSC \href{https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L363}{LDSC's munge_sumstats.py}
  # Reparation of B and SE from Z, P and EAF is from \href{https://www.biostars.org/p/319584/}
  cli::cli_h3("Repairing missing statistics columns:")
  cli::cli_ol()
  cli::cli_li("Transform OR to B if OR exists")
  cli::cli_li("Remove OR if both OR and B exists")
  cli::cli_li("Impute Z based on B and SE if both B and SE exist and Z is missing")
  cli::cli_li("Impute Z based on P and B if Z and SE is missing")
  cli::cli_li("Impute B and SE if both are missing, and Z, EAF and N is present")
  cli::cli_h3("Starting reparations:")
  start_cols <- colnames(tbl)


  if("OR" %in% colnames(tbl) & !"B" %in% colnames(tbl)) {
    cli::cli_alert_info("Found OR but not BETA. converting to B using {.code base::log(OR)}")
    tbl <- dplyr::mutate(tbl, B = log(OR)) |>
      dplyr::select(-OR)
  }

  if("OR" %in% colnames(tbl) & "B" %in% colnames(tbl)) {
    cli::cli_alert_info("found OR and B, removing OR")
    tbl <- dplyr::select(tbl, -OR)
  }

  if(all(c("B", "SE") %in% colnames(tbl)) & !"Z" %in% colnames(tbl)) {
    cli::cli_alert_info("Z missing: Calculating Z using the formula:  Z = B / SE")
    tbl <- dplyr::mutate(tbl, Z = B / SE)
  }

  if(all(c("B", "P") %in% colnames(tbl)) & !"Z" %in% colnames(tbl)) {
    cli::cli_alert_info("Found B and P but not Z. Imputing Z using:
                        sign(beta) * sqrt(stats::qchisq(pvalue,1,lower=FALSE))")
    tbl <- dplyr::mutate(tbl, Z = z_from_p_b(P, B))
  }

  if(all(c("Z", "N", "EAF") %in% colnames(tbl)) & all(!c("B", "SE") %in% colnames(tbl))) {
    cli::cli_alert_info("{.strong Imputing B and SE, using Z, EAF and N }. The B and SE will correspond to a standardized scale, which might not always be the same scale that the original betas was on.")
    # check variable N
    if(length(unique(tbl$N)) == 1) cli::cli_alert_danger("Attempting to repair B and SE from Z,EAF and N. However, N does not not vary across SNPs, indicating it's the study-wide N and not SNP-wise N. This will introduce issues with the conversion")
    tbl <- dplyr::mutate(tbl,
                         B = beta_from_z_eaf_n(Z,EAF,N),
                         SE = se_from_z_eaf_n(Z,EAF,N)
    )
  }
  if("Z" %in% colnames(tbl) & !"P" %in% colnames(tbl)){
    cli::cli_alert("Imputing P based on Z")
    tbl <- dplyr::mutate(tbl,P = stats::pnorm(-abs(tbl$Z)) *2)

  }



  # end message -------------------------------------------------------------

  end_cols <- colnames(tbl)
  new_cols <- end_cols[!end_cols %in% start_cols]
  if(length(new_cols) > 0) {
    cli::cli_h3("Finished repair_stats: ")
    cli::cli_alert_info("Added {length(new_cols)} new columns: {new_cols}")

  }


  # finished - return tibble ---------------------------------------------------
  tbl

}

flag_rsid_history <- function(tbl, rs_merge_arch) {



  # merge with RsMergeArch, and update where RSID matches old RSID
  dplyr::select(tbl, dplyr::all_of(c("rowid", "RSID"))) |>
    dplyr::left_join(rs_merge_arch, by = "RSID") |>
    # if RSID has been merged, update to correct one
    dplyr::mutate(new_RSID = dplyr::if_else(is.na(new_RSID), RSID, new_RSID)) |>
    dplyr::rename(old_RSID = RSID, RSID = new_RSID) |>
    dplyr::mutate(old_RSID = dplyr::if_else(RSID != old_RSID, old_RSID, NA_character_)) |>
    dplyr::select(dplyr::all_of(c("rowid", "old_RSID", "RSID")))


}

update_rsid <- function(tbl, merge_history) {
  stopifnot("RSID" %in% colnames(tbl))
  tmp <- dplyr::left_join(tbl, merge_history, by = "RSID") |>
    # update RSID where an updated RSID exists
    dplyr::mutate(new_RSID = dplyr::if_else(is.na(new_RSID), RSID, new_RSID)) |>
    # save old_RSID
    dplyr::rename(old_RSID = RSID, RSID = new_RSID) |>
    dplyr::mutate(old_RSID = dplyr::if_else(RSID != old_RSID, old_RSID, NA_character_))


}



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


beta_from_z_eaf_n <- function(Z, EAF, N) {
  # ref https://www.biostars.org/p/319584/
  # and https://www.nature.com/articles/ng.3538
  Z / sqrt((2*EAF*(1 - EAF)) * (N + (Z^2)))
}

se_from_z_eaf_n <- function(Z, EAF, N) {
  # ref https://www.biostars.org/p/319584/
  # https://www.nature.com/articles/ng.3538
  1/sqrt((2*EAF)*(1-(EAF))*(N+(Z^2)))
}

z_from_p_b <- function(pvalue, beta) {
  sign(beta) * sqrt(stats::qchisq(pvalue,1,lower=FALSE))
}

#' Fast detection of duplicates in RSID or CHR and POS columns
#' @description
#' creates new column: `dup_rsid` or `dup_chr_pos`, a T/F flag.
#' Specifically, flags both rows in a duplication pair, and not just first or
#' last duplicate row, making it easy to work with all rows that are part of a
#' duplication
#'
#' @param tbl a tibble with tidyGWAS column names
#' @param column construct id with either 'rsid' or 'chr_pos' or 'chr_pos_ref_alt' or 'rsid_ref_alt'
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' flag_duplicates(tbl, column = "chr_pos")
#' flag_duplicates(tbl, column = "crsid")
#'
#' }
flag_duplicates <- function(tbl, column = "rsid") {
  new_name <- glue::glue("dup_{column}")
  stopifnot("`column`  only takes 'rsid' or 'chr_pos' or 'chr_pos_ref_alt' or 'rsid_ref_alt' as arguments" = column %in% c("rsid", "chr_pos", "chr_pos_ref_alt", "rsid_ref_alt"))
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


