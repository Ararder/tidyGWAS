



guess_names <- function(tbl) {

  # find best mapping
  cols <- colnames(tbl)
  #
  cols <- cols[cols != "rowid"]
  scores <- purrr::map_dbl(formats, \(x) count_matches(cols, x))
  best_format_name <- names(scores)[which.max(scores)]
  best_format      <- formats[[best_format_name]]




  success <- best_format[best_format %in% cols]
  missing <- cols[!cols %in% best_format]

  cli::cli_h3("The detected format is {.val {best_format_name}}")
  cli::cli_inform("Was able to map {length(success)} out of {length(cols)} to the {.val {best_format_name}} format")
  for(i in seq_along(success)) cli::cli_alert_success("{names(success)[i]} -> {success[i]}")

  if(length(missing) > 0) cli::cli_alert_warning("Failed to map: {.val {missing}}")

  if(length(missing)) {
    cli::cli_alert_info("Attempting to map: {.val {missing}} to a dictionary of column names")

    mapped <- purrr::map(missing,\(x) find_match(x, colname_dict))

    has_entry <- purrr::map_lgl(mapped, \(x) !rlang::is_empty(x))


    succesful_mapping <- missing[has_entry] |>
      purrr::set_names(mapped[has_entry])

    if(length(succesful_mapping) > 0) {
      for(i in seq_along(succesful_mapping)) cli::cli_alert_success("{names(succesful_mapping)[i]} -> {succesful_mapping[i]}")
      success <- c(success, succesful_mapping)
    } else {
      cli::cli_alert_danger("No additional mappings found")
    }

  }




  missing <- setdiff(cols, unname(unlist(success)))
  # rowid should not be reported as an unmapped column
  missing <- missing[missing != "rowid"]
  if (length(missing)) cli::cli_alert_warning(
    "extra column(s) {.val {missing}} were not mapped to a tidyGWAS column"
  )



  # finished ----------------------------------------------------------------
  dplyr::rename(tbl, !!!success)



}

count_matches <- function(cnames, fmt) {
  fmt_cols <- unname(unlist(fmt))
  sum(fmt_cols %in% cnames)
}

find_match <- function(name, dict) {

  names(dict)[purrr::map_lgl(dict, \(potential_col_names) name %in%  potential_col_names)]

}


