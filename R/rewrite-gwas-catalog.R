
#' Download summary statistics from GWAS catalog
#'
#' @param study_id A single character string with study ID, e.g. "GCST90475332"
#' @param quiet TRUE/FALSE - controls progress bar for downloading
#'
#' @returns a filepath to downloaded summary statistics
#' @export
#'
#' @examples \dontrun{
#' from_gwas_catalog("GCST90475332")
#' }
from_gwas_catalog <- function(study_id,quiet = FALSE) {



  rlang::check_required(study_id)
  rlang::is_scalar_logical(quiet) ||
    cli::cli_abort("quiet: {.arg {quiet}} must be a single logical value")
  rlang::is_scalar_character(study_id) ||
    cli::cli_abort(
      "study_id: {.arg {study_id}} must be a single character string"
    )
  workdir <- fs::dir_create(fs::path(tempdir(), study_id))

  # -------------------------------------------------------------------------
  # 2. Get URLs
  ftp_base <- get_ftp_path(study_id)
  cli::cli_alert_info("Checking: {.url {ftp_base}}")
  target_file <- pick_file(ftp_base, study_id = study_id)

  # -------------------------------------------------------------------------
  # 3. Download
  dest_file <- fs::path(workdir, basename(target_file))
  curl::curl_download(target_file, dest_file, quiet = quiet)

  dest_file


}

pick_file <- function(ftp_base, study_id) {

  target_file <- tryCatch({
    # Try Harmonised first
    h_url <- paste0(ftp_base, "harmonised/")
    files <- .list_ftp_files(h_url)

    # Filter for standard harmonised extensions
    valid <- files[grepl("\\.h\\.tsv\\.gz$|\\.f\\.tsv\\.gz$", files)]

    if (length(valid) > 0) {
      cli::cli_alert_success("Found Harmonised data")
      paste0(h_url, valid[1])
    } else {
      stop("No harmonised files")
    }
  }, error = function(e) {
    # Fallback to Raw directory
    cli::cli_alert_warning("Harmonised not found. Checking raw directory...")
    files <- .list_ftp_files(ftp_base)

    # Filter: ends in tsv/txt/gz, ignore metadata/images
    valid <- files[grepl("\\.(tsv|txt|tsv\\.gz|txt\\.gz)$", files, ignore.case = TRUE)]
    valid <- valid[!grepl("readme|meta|md5", valid, ignore.case = TRUE)]

    if (length(valid) == 0) cli::cli_abort("No valid summary statistics found for {study_id}")

    # Heuristic: Pick the one containing the ID, or the first one
    best <- valid[grepl(study_id, valid)]
    if (length(best) == 0) best <- valid[1]

    cli::cli_alert_info("Selected raw file: {best}")
    paste0(ftp_base, best)
  })

  target_file
}

#' Internal helper to list FTP files cleanly
.list_ftp_files <- function(url) {
  h <- curl::new_handle(dirlistonly = TRUE)
  con <- try(curl::curl_fetch_memory(url, handle = h), silent = TRUE)

  if (inherits(con, "try-error") || con$status_code >= 400) {
    cli::cli_abort("Could not find a directory corresponding to the provided GCST ID.
                   Are you sure summary statistics are available for the GCST ID?")
  }

  # FTP listing comes as a newline separated string
  content <- rawToChar(con$content)
  files <- strsplit(content, "[\r\n]+")[[1]]
  return(files)
}

get_ftp_path <- function(study_id) {
  num_part <- as.numeric(sub("GCST", "", study_id))

  # Calculate range start
  range_start <- (floor((num_part - 1) / 1000) * 1000) + 1
  range_end   <- range_start + 999

  # Determine width based on the input string length (robustness fix)
  # If input is GCST90001234 (12 chars), digits are 8.
  width <- nchar(study_id) - 4
  fmt   <- paste0("GCST%0", width, "d")

  range_str <- paste0(sprintf(fmt, range_start), "-", sprintf(fmt, range_end))
  paste0("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/", range_str, "/", study_id, "/")
}

