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
from_gwas_catalog <- function(study_id, quiet = FALSE, harmonised = FALSE) {
  rlang::check_required(study_id)
  rlang::is_scalar_logical(harmonised) || cli::cli_abort("harmonised: {.arg {harmonised}} must be a single logical value")
  rlang::is_scalar_logical(quiet) || cli::cli_abort("quiet: {.arg {quiet}} must be a single logical value")
  rlang::is_scalar_character(study_id) || cli::cli_abort("study_id: {.arg {study_id}} must be a single character string")
  check_sumstats_avail(study_id) ||cli::cli_abort("Full summary statistics for {study_id} must be available in GWAS catalog")


  workdir <- fs::dir_create(fs::path(tempdir(), study_id))

  all_urls <- scrape_dir(get_study_id_url(study_id))
  picked <- .pick_sumstats(all_urls, harmonised)


  meta_file <- fs::path(workdir, fs::path_file(picked$yaml))
  gwas_file <- fs::path(workdir, fs::path_file(picked$gwas))

  # download meta file and print it?


  cli::cli_alert_info("Downloading GWAS summary statistics {.path {gwas_file}}")
  curl::curl_download(url = picked$yaml, destfile = meta_file)
  curl::curl_download(url = picked$gwas, destfile = gwas_file, quiet = quiet)

  gwas_file
}


get_study_id_url <- function(study_id) {
  root <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

  num <- as.integer(sub("GCST", "", study_id))
  lo <- ((num - 1) %/% 1000) * 1000 + 1
  hi <- lo + 999

  range_folder <- sprintf("GCST%08d-GCST%08d", lo, hi)
  file.path(root, range_folder, study_id)
}


scrape_dir <- function(url, pattern = "\\.(gz|tsv|yaml|tbi|log|txt)$") {
  page <- xml2::read_html(url)
  hrefs <- rvest::html_attr(rvest::html_nodes(page, "a"), "href")
  hrefs <- hrefs[!grepl("^\\?|^\\.{2}/|^/pub|^$", hrefs)]
  non_harmonised <- file.path(url, hrefs[!grepl("harmonised/", hrefs)])

  if (any(grepl("harmonised/", hrefs))) {
    page <- xml2::read_html(file.path(url, "harmonised/"))
    hrefs <- rvest::html_attr(rvest::html_nodes(page, "a"), "href")
    hrefs <- hrefs[!grepl("^\\?|^\\.{2}/|^/pub|^$", hrefs)]

    c(non_harmonised, file.path(file.path(url, "harmonised"), hrefs))
  } else {
    non_harmonised
  }
}


.pick_sumstats <- function(urls, harmonised = TRUE) {
  ## 1. keep only real data files (drop .yaml, .tbi, .log, etc.)


  if (!length(urls)) {
    cli::cli_abort(
      "No .tsv/.txt summary-statistics files found in the directory"
    )
  }

  harmonised_files <- fs::path_dir(urls) |> fs::path_file() == "harmonised"
  if (!any(harmonised_files) & harmonised) {
    cli::cli_abort(
      "No harmonised summary-statistics files found in the directory"
    )
  } else if(harmonised) {

    gwas <- stringr::str_subset(urls[harmonised_files], "\\.tsv(\\.gz)?$")
    yaml <- stringr::str_subset(urls[harmonised_files], "meta.yaml$")
  } else if(!harmonised) {
    gwas <- stringr::str_subset(urls[!harmonised_files], "\\.tsv(\\.gz)?$")
    yaml <- stringr::str_subset(urls[!harmonised_files], "meta.yaml$")

  }

  list(
    gwas = gwas,
    yaml = yaml
  )


}


# -------------------------------------------------------------------------

get_gwas_catalog_region <- function(study_id, chr, start, end) {
  rlang::is_scalar_character(study_id) ||
    cli::cli_abort(
      "study_id: {.arg {study_id}} must be a single character string"
    )
  endpoint <- paste0(base_url, "/chromosomes/", chr, "/associations")
  base_url <- "https://www.ebi.ac.uk/gwas/summary-statistics/api"

  params <- list(
    bp_lower = start,
    bp_upper = end,
    study_accession = study_id,
    size = 5000
  )

  t <- .fetch_page(endpoint, params)

  tidy_catalog_json(t)
}


check_sumstats_avail <- function(
  study_id,
  server = "https://www.ebi.ac.uk/gwas/rest/api"
) {
  resp <- httr::GET(
    paste0(server, "/studies/", study_id),
    httr::accept("application/json")
  )
  httr::stop_for_status(resp)
  jsonlite::fromJSON(httr::content(
    resp,
    "text",
    encoding = "UTF-8"
  ))$fullPvalueSet
}


check_rest_avail <- function(study_id) {
  url <- sprintf(
    "https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/%s",
    study_id
  )
  httr::status_code(httr::GET(url, httr::user_agent("has_rest_api"))) == 200
}


# -------------------------------------------------------------------------

.fetch_page <- function(url, query) {
  resp <- httr::GET(url, query = query)
  httr::stop_for_status(resp)

  data <- jsonlite::fromJSON(
    httr::content(resp, "text", encoding = "UTF-8"),
    flatten = TRUE
  )
}


variables <- c(
  "variant_id",
  "base_pair_location",
  "effect_allele",
  "other_allele",
  "beta",
  "se",
  "effect_allele_frequency",
  "odds_ratio",
  "ci_lower",
  "ci_upper",
  "p_value"
)

tidy_catalog_json <- function(data) {
  assocs <- data[["_embedded"]][["associations"]]
  purrr::map(assocs, \(row) {
    existing <- names(row)
    does_exist <- intersect(variables, existing)

    row[does_exist] |>
      purrr::map(\(x) ifelse(is.null(x), NA, x)) |>
      dplyr::as_tibble()
  }) |>
    purrr::list_rbind()
}
