utils::globalVariables(c("variant.chromosome", "variant.position",
                        "variant.referenceAllele", "variant.alternateAllele",
                        "study.traitFromSource", "study.id", "study.diseases",
                        "locus.rows", "l2GPredictions.rows", "locusSize.count",
                        "pValueMantissa", "pValueExponent"))
remove_cols <- c(
  "studyLocusId",
  "finemappingMethod",
  "confidence",
  "study.diseases",
  "locus.rows",
  "l2GPredictions.rows",
  "locusSize.count",
  "variant.id",
  "pValueMantissa",
  "pValueExponent"
)


#' Query Open targets for all credible sets containing the variant
#'
#' @param variant_id in the open targets format chr_pos_ref_alt: 5_100_101_A_C
#' @param page_size number of results per page, default 500
#' @param api_url URL of the Open Targets GraphQL API, default
#'
#' @returns a [dplyr::tibble] with the following columns:
#'  - `CHR`: chromosome
#'  - `POS`: position
#'  - `P`: p-value
#'  - `REF`: reference allele
#'  - `ALT`: alternate allele
#'  - `description`: trait description
#'  - `study_id`: study ID
#'
#' @export
#'
#' @examples \dontrun{
#' get_open_targets_cs("7_140459051_C_G")
#' }
get_open_targets_cs <- function(
  variant_id,
  page_size = 500,
  api_url = "https://api.platform.opentargets.org/api/v4/graphql"
) {
  rlang::is_scalar_character(variant_id) ||
    cli::cli_abort(
      "Incorrect value for `variant_id`: {variant_id}.
                   Must be a single character string."
    )

  # ---- first request (gets `count`) -----------------------------------------

  pg0 <- fetch_page(
    variant_id = variant_id,
    page_size = page_size,
    page_index = 0,
    api_url = api_url
  )
  info <- pg0$variant$gwasCredibleSets
  total <- info$count %||% 0L

  # handle pagination
  rows <- info$rows
  if (total > page_size) {
    extra_pages <- seq(page_size, total - 1L, by = page_size) / page_size
    for (i in extra_pages) {
      pg <- fetch_page(i)
      rows <- vctrs::vec_c(rows, pg$variant$gwasCredibleSets$rows)
    }
  }

  # ---- tidy up --------------------------------------------------------------
  rows |>
    dplyr::as_tibble() |>
    dplyr::mutate(P = pValueMantissa * 10^pValueExponent) |>
    dplyr::select(
      -dplyr::any_of(remove_cols)
    ) |>
    dplyr::relocate(
      CHR = variant.chromosome,
      POS = variant.position,
      P,
      REF = variant.referenceAllele,
      ALT = variant.alternateAllele,
      description = study.traitFromSource,
      study_id = study.id
    )
}


fetch_page <- function(variant_id, page_size, page_index = 0,api_url) {
  # construct query body
  body <- list(
    query = gql_variant_template,
    variables = list(
      variantId = variant_id,
      size = page_size,
      index = page_index
    )
  )
  res <- httr::POST(
    api_url,
    body = body,
    encode = "json",
    httr::user_agent("otg-graphql-r/0.1")
  )
  httr::stop_for_status(res)

  jsonlite::fromJSON(
    httr::content(res, "text", encoding = "UTF-8"),
    flatten = TRUE
  )[["data"]]
}


# ---- GraphQL query --------------------------------------------------------
gql_variant_template <- '
    query GWASCredibleSetsQuery($variantId: String!, $size: Int!, $index: Int!) {
      variant(variantId: $variantId) {
        id
        referenceAllele
        alternateAllele
        gwasCredibleSets: credibleSets(studyTypes: [gwas],
                                       page: { size: $size, index: $index }) {
          count
          rows {
            studyLocusId
            pValueMantissa
            pValueExponent
            beta
            finemappingMethod
            confidence
            variant {
              id
              chromosome
              position
              referenceAllele
              alternateAllele
            }
            study {
              traitFromSource
              id
              diseases {
                name
                id
                therapeuticAreas {
                  name
                  id
                }
              }
            }
            locus(variantIds: [$variantId]) {
              rows { posteriorProbability }
            }
            locusSize: locus { count }
            l2GPredictions(page: { size: 1, index: 0 }) {
              rows {
                target   { id approvedSymbol }
                score
              }
            }
          }
        }
      }
    }'
