# -------------------------------------------------------------------------

#' Query a specific region of interest for a using a gwas catalog study_id
#'
#' @param study_id a study accession ID, e.g. "GCST000001"
#' @param chr chromosome number, e.g. "1" - not "chr1"
#' @param start base pair start position, e.g. 1000000
#' @param end base pair end position, e.g. 2000000
#'
#' @returns a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' #' get_gwas_catalog_region("GCST000001", "1", 1000000, 2000000)
#' }
from_gwas_catalog_region <- function(study_id, chr, start, end) {
  rlang::is_scalar_character(study_id) ||
    cli::cli_abort(
      "study_id: {.arg {study_id}} must be a single character string"
    )

  check_rest_avail(study_id) ||
    cli::cli_abort(
      "REST API is not available for study {.val {study_id}}.
      The study might not have been `ingested` yet, or summary statistics are not available.
      ")

  base_url <- "https://www.ebi.ac.uk/gwas/summary-statistics/api"
  endpoint <- paste0(base_url, "/chromosomes/", chr, "/associations")

  params <- list(
    bp_lower = start,
    bp_upper = end,
    study_accession = study_id,
    size = 5000
  )

  t <- .fetch_page(endpoint, params)

  tidy_catalog_json(t)
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
