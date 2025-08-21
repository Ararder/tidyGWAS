utils::globalVariables(c("ancestry", "est"))
#' Estimate ancestry composition from allele frequencies
#'
#'
#' This function implements the `bigsnpr::snp_ancestry_summary()`
#' written by Florian Prive.
#'
#' The function and reference data is based on
#' <https://doi.org/10.1093/bioinformatics/btac348>
#'
#' [ancestry_comp()] requires additional reference data. This data is
#' only available if you have downloaded the updated reference tidyGWAS data
#' The reference file `ancestry_data.parquet` is expected to exist in the
#' same dbSNP folder that is used with [tidyGWAS()]
#'
#' @param tbl a data.frame with atleast columns `RSID`, `EffectAllele`, `OtherAllele`, and `EAF`.
#' @param dbsnp_path filepath to the dbSNP reference data. Same path as for [tidyGWAS()]
#' @param min_cor minimum correlation between predicted and observed allele frequencies.
#' @param sum_to_one Force ancestry coefficients to sum to 1?
#'
#' @returns a tibble with columns `est` and `ancestry`
#' @export
#'
#' @examples \dontrun{
#' ancestry_comp(tbl, "path_to_data")
#'
#' }
ancestry_comp <- function(tbl, dbsnp_path, min_cor = 0.4, sum_to_one=TRUE) {
  # NOTE: this is almost a straight copy of the function from the bigsnpr package,
  # all credit to @privefl
  rlang::check_required(tbl)
  rlang::is_scalar_double(min_cor) || cli::cli_abort("`min_cor` must be a single numeric value.")
  rlang::is_scalar_logical(sum_to_one) || cli::cli_abort("`sum_to_one` must be a single logical value.")
  rlang::check_installed("quadprog")
  rlang::check_installed("Matrix")

  ref_data_path <- file.path(dbsnp_path, "ancestry_data.parquet")
  if(!file.exists(ref_data_path)) {
    cli::cli_abort("Reference data not found at {ref_data_path}. Please download it using `get_data(dbsnp_path)`")
  }

  # -------------------------------------------------------------------------
  data <- format_anc_comp_data(tbl, path = ref_data_path)

  projection <- data$projections
  correction <- data$correction
  freq <- data$gwas_freq
  X0 <- data$ref_freq
  if (mean(stats::cor(X0, freq)) < -0.2) {
    cli::cli_abort("Frequencies seem all reversed; switch reference allele?")
  }
  # project allele frequencies onto the PCA space
  X <- crossprod(projection, X0)
  y <- crossprod(projection, freq) * correction

  cp_X_pd <- Matrix::nearPD(crossprod(X), base.matrix = TRUE)
  if (!isTRUE(cp_X_pd$converged)) {
    cli::cli_abort("Could not find nearest positive definite matrix.")
  }

  # solve QP problem using https://stats.stackexchange.com/a/21566/135793
  res <- quadprog::solve.QP(
    Dmat = cp_X_pd$mat,
    dvec = crossprod(y, X),
    Amat = cbind(-1, diag(ncol(X))),
    bvec = c(-1, rep(0, ncol(X))),
    meq  = if (sum_to_one) 1 else 0
  )

  cor_pred <- drop(stats::cor(drop(X0 %*% res$solution), freq))
  if (cor_pred < min_cor) {
    cli::cli_abort("Correlation between frequencies is too low: {cor_pred}, check matching between variants.")
  }
  if (cor_pred < 0.99) {
    cli::cli_alert_warning("The solution does not perfectly match the frequencies.")

  }



  dplyr::tibble(
    est = round(res$solution, 7),
    ancestry = colnames(X0),
    allele_freq_corr = drop(stats::cor(X0, freq)),
    cor_pred = cor_pred
  ) |>
    # dplyr::arrange(dplyr::desc(est))
    dplyr::mutate(
      # recode ancestry as in tutorial: https://privefl.github.io/bigsnpr/articles/ancestry.html
      ancestry = dplyr::case_when(
        ancestry %in% c("scandinavia", "united_kingdom", "ireland") ~ "europe_north_west",
        ancestry %in% c("europe_south_east", "europe_north_east") ~ "europe_east",
        .default = ancestry
      )
    ) |>
    dplyr::arrange(dplyr::desc(est)) |>
    dplyr::summarise(est = sum(est), .by = "ancestry")


}




# -------------------------------------------------------------------------


format_anc_comp_data <- function(tbl, path) {
  ref_data <- arrow::read_parquet(path)
  pcs_names <- paste0("pc", seq(1:16))


  correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                  1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

  freq_ref_names <- colnames(ref_data)[4:24]


  # -------------------------------------------------------------------------
  # first filter, remove NA's

  snps <- dplyr::select(tbl, dplyr::all_of(c("RSID", "EffectAllele", "OtherAllele", "EAF"))) |>
    tidyr::drop_na()
  cli::cli_alert_info("dropped {nrow(tbl) - nrow(snps)} rows with missing values.")



  # -------------------------------------------------------------------------
  # second filter, intersection of variants

  one <- dplyr::inner_join(snps, ref_data, by = c("RSID", "EffectAllele", "OtherAllele")) |> dplyr::collect()
  two <- dplyr::inner_join(
    snps,
    ref_data,
    by = c("RSID", "EffectAllele" = "OtherAllele", "OtherAllele" = "EffectAllele")
    ) |>
    dplyr::mutate(EAF = 1-EAF)

  merged <- dplyr::bind_rows(one,two)

  cli::cli_alert_success(
    "Started with {nrow(snps)} variants in GWAS data and {nrow(ref_data)} in reference data, found {nrow(merged)} matching variants."
    )
  cli::cli_alert_info(
    "{nrow(two)} variants were reversed to match effect allele frequency, {nrow(one)} were kept as is."
  )

  # return format
  list(
    gwas_freq = merged$EAF,
    ref_freq = dplyr::select(merged, dplyr::all_of(freq_ref_names)) |> as.matrix(),
    projections = dplyr::select(merged, dplyr::all_of(pcs_names)) |> as.matrix(),
    correction = correction
  )

}

