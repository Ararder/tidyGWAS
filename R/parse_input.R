parse_input <- function(tbl, ...) {
  # -------------------------------------------------------------------------

  prevent_err(...)
  args <- rlang::list2(...)
  if("delim" %in% names(args)) {
    cli::cli_alert_warning(
      "
          The use of the delim argument has been superceded.
          tidyGWAS now uses data.table::fread() for parsing files.
          In most instances this should result in automatic detection of the delimiter.
          Arguments can still be passed from tidyGWAS() to fread() with ..."
    )
  }



  # -------------------------------------------------------------------------
  if(is_data_frame(tbl)) {
    gwas <- tbl
  } else if (is_gwas_catalog(tbl)) {

    cli::cli_alert_success(
      "Detected a study ID: {.emph {tbl}}. Fetching data from the GWAS Catalog"
    )
    fp <- from_gwas_catalog(study_id = tbl)
    gwas <- arrow::read_delim_arrow(fp, delim = "\t")

  } else if(is_url(tbl)) {


    cli::cli_alert_info(
      "Detected a URL: {.emph {tbl}}. Attempting to download the file."
    )
    tempdir <- tempdir()
    dest <- file.path(tempdir, basename(tbl))


    curl::curl_download(
      url = tbl,
      destfile = dest,
      quiet = FALSE,
    )

    gwas <- data.table::fread(dest)


  } else if (is_local_filepath(tbl)) {

    gwas <- data.table::fread(tbl, ...)
  }




  if(rlang::is_scalar_character(tbl)) {
    filename <- basename(tbl)
    md5 <- tools::md5sum(tbl)
  } else {
    filename <- "raw"
    md5 <- NULL
  }


  if (!"rowid" %in% colnames(gwas)) {
    gwas$rowid <- 1:nrow(gwas)
  }

  list("tbl" = gwas, "filename" = filename, "md5" = md5)
}



#
is_gwas_catalog <- function(char) {
  stringr::str_detect(char, "^GCST\\d+$")
}



#
is_url <- function(x) {


  url_regex <- paste0(
    "^(?:[A-Za-z][A-Za-z0-9+.-]{1,}://)",
    "(?:[^@/\\s]+@)?",
    "(?:\\[[0-9a-fA-F:]+\\]|",
    "(?:[A-Za-z0-9.-]+))",
    "(?::[0-9]{1,5})?",
    "(?:[/?#][^\\s]*)?$"
  )

  result <- grepl(url_regex, x, perl = TRUE)

  # Treat NA or empty strings as FALSE
  result[is.na(x) | x == ""] <- FALSE
  result

}


#
is_data_frame <- function(x) {
  "data.frame" %in% class(x)

}



#
is_local_filepath <- function(obj) {
  file.exists(obj)

}

prevent_err <- function(delim) {
  return(NULL)
}

