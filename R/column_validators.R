utils::globalVariables(c("P", "p_was_0", "B", "EAF", "N", "missing_ea_oa", "SE", "non_acgt", "Z", "indel", "OR"))
validate_rsid <- function(tbl) {


  # setup -------------------------------------------------------------------

  start_message("RSID")
  tbl$RSID <- as.character(tbl$RSID)
  tbl <- flag_incorrect_rsid_format(tbl)
  invalid_rsid_format <- dplyr::filter(tbl, invalid_rsid)


  # No invalid RSIDs --------------------------------------------------------

  if(nrow(invalid_rsid_format) == 0) {
    cli::cli_alert_info("Found no rows with invalid RSID format")
    return(list("data" = tbl,"chr_pos" =  NULL,"failed" = NULL))
  }


  # invalid RSIDs, but other columns exist ----------------------------------

  if(all(c("CHR", "POS", "EffectAllele", "OtherAllele") %in% colnames(tbl))) {
    cli::cli_alert_info("Found { nrow(dplyr::filter(tbl, invalid_rsid)) } rows with invalid RSID format. ")
    return(
      list(
        "data" = tbl,
        "chr_pos" = dplyr::filter(tbl, invalid_rsid) |> dplyr::select(-RSID, -invalid_rsid),
        "failed" = NULL)
      )
  }


  # no CHR or POS, so need to parse RSID ------------------------------------

  cli::cli_alert_info("Found {nrow(invalid_rsid_format)} rows with invalid RSID format: ")
  cli::cli_alert_info("Attempting to parse format...")


  # check if its possible to get CHR:POS:REF:ALT or CHR:POS from invalid RSIDs
  # UPDATE: i think it's much better to always use EffectAllele and OtherAllele
  # from the sumstats, than trying to infer it, as it's unknown which is the
  # effect allele when parsing CHR:POS:REF:ALT
  attempt <- split_rsid_by_regex(invalid_rsid_format) |>
    dplyr::select(-EffectAllele, -OtherAllele) |>
    dplyr::inner_join(dplyr::select(tbl, rowid, EffectAllele, OtherAllele), by = "rowid")


  # print results to user, if any rows have been parsed correctly
  if(nrow(attempt) > 0) cli::cli(c(cli::cli_alert_info("Parsed format"), cli::cat_print(head(dplyr::select(chr_pos, RSID, CHR, POS, EffectAllele, OtherAllele), 5), file = stderr())))

  # check if failed to parse any rows
  failed <- dplyr::anti_join(dplyr::filter(tbl, invalid_rsid), attempt, by = "rowid")

  list("data" = tbl, "chr_pos" =  attempt,"failed" = failed)

}


validate_chr <- function(tbl) {

  start_message("CHR")
  valid_chr <- c(1:22, "X", "Y", "MT")
  tbl <- dplyr::mutate(tbl,
    CHR = as.character(CHR),
    CHR = stringr::str_to_upper(CHR),
    # can sometimes be chr22, or ch22
    CHR = stringr::str_remove(CHR, "CHR"),
    CHR = stringr::str_remove(CHR, "CH"),
    # can now handle
    CHR = dplyr::if_else(CHR == "23", "X", CHR),
    CHR = dplyr::if_else(CHR == "M", "MT", CHR),
    invalid_chr = dplyr::if_else(!CHR %in% valid_chr | is.na(CHR), TRUE, FALSE)
    )

  end_message(tbl, "chr")
  tbl

}





validate_pos <- function(tbl) {
  start_message("POS")

  tbl <- dplyr::mutate(tbl,
      POS = as.integer(POS),
      invalid_pos = dplyr::if_else(POS <= 0 | !is.finite(POS) | POS >= 10^9, TRUE, FALSE)
  )
  end_message(tbl, "pos")
  tbl


}

validate_P <- function(tbl, convert_p = 2.225074e-308) {
  start_message("P", convert_p = convert_p)


  tbl <-
    dplyr::mutate(tbl,
      P = as.double(P),
      p_was_0 = dplyr::if_else(P == 0, "Yes", "No"),
      P = dplyr::if_else(p_was_0 == "Yes", convert_p, P),
      invalid_P = dplyr::if_else(!is.finite(P) | P > 1 | P < 0, TRUE, FALSE),
      ) |>
    dplyr::select(-p_was_0)

  end_message(tbl, "P")


  tbl

}


validate_ea_oa <- function(tbl) {
  stopifnot(all(c("EffectAllele", "OtherAllele") %in% colnames(tbl)))
  stopifnot(
    "EffectAllele and OtherAllele are not character columns, indicating mislabelling of columns" =
      is.character("EffectAllele") & is.character("OtherAllele")
  )
  start_message("EA_OA")
  # cli::cli_alert_info("Checking EffectAllele and OtherAllele: Only values in [A,C,G,T] are allowed")
  # cli::cli_alert_info("Converting EffectAllele and OtherAllele to uppercase")
  cols_to_start <- colnames(tbl)

  tmp <- tbl |>
    dplyr::mutate(
      EffectAllele = stringr::str_to_upper(EffectAllele),
      OtherAllele = stringr::str_to_upper(OtherAllele)
    )
  possible_alleles <- c("A", "C","G","T")
  ea_alleles <- dplyr::count(tmp, EffectAllele)
  oa_alleles <- dplyr::count(tmp, OtherAllele)

  if(any(!ea_alleles$EffectAllele %in% possible_alleles)) {
    cli::cli_alert_danger("Found alleles other than A,C,G,T in EffectAllele: ")
    cli::cat_print(ea_alleles)

  }

  if(any(!oa_alleles$OtherAllele %in% possible_alleles)) {
    cli::cli_alert_danger("Found alleles other than A,C,G,T in OtherAllele: ")
    cli::cat_print(oa_alleles)
  }


  tmp <- dplyr::mutate(tmp, invalid_ea_oa = dplyr::if_else(!EffectAllele %in% possible_alleles | !OtherAllele %in% possible_alleles, TRUE, FALSE))
  n_invalid <- sum(tmp$invalid_ea_oa)

  if(n_invalid > 0) {
    cli::cli_alert_warning("Found {n_invalid} rows with non-ACGT codes")
  } else {
    cli::cli_alert_success("All rows pass EA/OA validation")
  }

  tmp

}



validate_b <- function(tbl) {
  start_message("B")
  tbl <- dplyr::mutate(tbl, B = as.double(B))


  # check median value ------------------------------------------------------

  median <- median(tbl$B, na.rm = TRUE)
  if(dplyr::between(median, left = 0.9, right = 1.1)) {
    cli::cli_alert_danger("WARNING: The median value of B is {median}, indicating that B has been mislabelled and contains odds-ratios")
  } else if(abs(median) > 0.1) {
    cli::cli_alert_danger("WARNING: The median value of B is {median}, which seems high")
  } else {
    cli::cli_alert_info("The median value of B is {median}, which seems reasonable")
  }


  # check for non finite ----------------------------------------------------

  tbl <- dplyr::mutate(tbl, invalid_B = dplyr::if_else(!is.finite(B), TRUE, FALSE))

 end_message(tbl, "B")

  tbl
}

validate_se <- function(tbl) {
  start_message("SE")

  tbl <- dplyr::mutate(tbl,
    SE = as.double(SE),
    invalid_SE = dplyr::if_else(SE <= 0, TRUE, FALSE)
    )

  end_message(tbl, "SE")
  tbl

}


validate_eaf <- function(tbl) {
  stopifnot("EAF" %in% colnames(tbl))
  start_message("EAF")


  tbl <-dplyr::mutate(tbl, EAF = as.double(EAF), invalid_EAF = dplyr::if_else(EAF <= 0 | EAF >= 1 | !is.finite(EAF), TRUE, FALSE))
  end_message(tbl, "EAF")

  tbl

}

validate_n <- function(tbl) {
  stopifnot("N" %in% colnames(tbl))
  start_message("N")

  tbl <- dplyr::mutate(tbl,N = as.integer(N),invalid_N = dplyr::if_else(N <= 0 | !is.finite(N), TRUE, FALSE))


  end_message(tbl, "N")
  tbl

}

validate_z <- function(tbl) {
  stopifnot("Z" %in% colnames(tbl))
  tbl <- dplyr::mutate(tbl, Z = as.double(Z))
  maxval <- max(abs(tbl$Z))

  if(maxval < 100) {
    cli::cli_alert_info("Found {maxval} as largest absolute Z score, which seems reasonable")
  } else {
    cli::cli_alert_warning("WARNING: Found {maxval} as largest absolute Z score, which seems highy unlikely")
  }


  tbl <- dplyr::mutate(tbl, invalid_Z = dplyr::if_else(!is.finite(Z), TRUE, FALSE))

  end_message(tbl, "Z")


  tbl
}


start_message <- function(col, convert_p) {
  if(col == "P") {

    cli::cli_h3("Validating the P column:")
    cli::cli_ol()
    cli::cli_li("Coerce to double {.code base::as.double()}")
    cli::cli_li("If P == 0, P is converted to {.arg {convert_p}}. Adjust by passing convert_p = 'your_preferred_value'")
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

  if(col == "EA_OA") {
    cli::cli_h3("Validating EffectAllele and OtherAllele:")
    cli::cli_ol()
    cli::cli_li("Will error if not type = character")
    cli::cli_li("EffectAllele and OtherAllele converted to uppercase")
    cli::cli_li("Check for EA/OA values that are not A,C,G or T")
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

  n_invalid <- sum(tbl[[glue::glue("invalid_{col}")]])
  if(n_invalid > 0) {
    cli::cli_alert_warning("{n_invalid} rows failed {col} validation")
  } else {
    cli::cli_alert_success("All rows pass {col} validation")

  }

}
