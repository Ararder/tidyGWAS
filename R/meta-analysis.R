utils::globalVariables(
  c("REF", "EA_is_ref", "CHR", "ID", "W", "N_info", "INFO","tmp", "EffectiveN", "N_EAF")
)


#' Perform meta-analysis of GWAS summary statistics datasets cleaned by tidyGWAS
#'
#' @param dset an [arrow::open_dataset()] object
#' @param method method to use for performing meta-analysis. Currently, only IVW (based on standard errors) is supported.
#' @param by a character vector of column names to group by. Default is c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele")
#'  Columns not in by will not be kept. So if you group by c("CHR", "POS"), only these columns will be kept in the output.
#'  The columns passed are used to define a unique variant. If for example you pass c("CHR", "POS"), then all variants with the same
#'  chromosome and position would be considered the same variant. You can use 'ID'
#'  to speed up the meta analysis, but RSID information will then be lost.
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' dset <- arrow::open_dataset("path_to/sumstats/")
#' res <- meta_analyze(dset)
#' }
#'
meta_analyze <- function(dset, by = c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele"), method = c("ivw")) {
  method <- rlang::arg_match(method)
  purrr::map(c(1:22), \(chrom) meta_analyze_by_chrom(dset, chrom = chrom, by = by)) |>
    purrr::list_rbind()
}


align_to_ref <- function(dset) {
  # EffectAllele is harmonized to always be the reference allele
  dset |>
    dplyr::mutate(
      EA_is_ref = dplyr::if_else(EffectAllele == REF, TRUE,FALSE),
      tmp = EffectAllele,
      EffectAllele = dplyr::if_else(EA_is_ref, EffectAllele, OtherAllele),
      OtherAllele = dplyr::if_else(EA_is_ref, OtherAllele, tmp),
      B = dplyr::if_else(EA_is_ref, B, B*-1),
      EAF = dplyr::if_else(EA_is_ref, EAF, 1-EAF)
    ) |>
    dplyr::select(-dplyr::all_of(c("EA_is_ref", "tmp")))

}

#' meta_analyze summary statistics, one chromosome at a time!
#' This function is exposed to allow for testing using real data
#' @inheritParams meta_analyze
#' @param chrom chromosome to use for meta-analysis
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' meta_analyze_by_crom(dset, chrom = "22")
#' }
#'
meta_analyze_by_chrom <- function(dset, chrom, by) {
  stats <- c("B", "SE", "EAF", "N", "CaseN", "ControlN","INFO", "EffectiveN")
  cols <- c(by, stats)

  dset |>
    dplyr::filter(CHR == {{ chrom }}) |>
    dplyr::filter(is.na(indel) | !indel) |>
    # dplyr::filter(is.finite(B) & is.finite(SE)) |>
    align_to_ref() |>
    dplyr::select(dplyr::any_of(cols)) |>
    # for each varianat, calculate the weight, and multiply B by weight
    # allele frequency and INFO is weighted by sample size
    dplyr::mutate(
      W = 1 / (SE^2),
      B = B*W,
      dplyr::across(dplyr::any_of(c("EAF", "INFO")), ~.x * N)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_contributions = dplyr::n(),
      dplyr::across(dplyr::any_of(c("W", "B")), sum),
      dplyr::across(dplyr::any_of(c("EAF", "INFO", "CaseN", "ControlN", "N", "EffectiveN")), ~sum(.x, na.rm=T)),
      # calculate the total sample size for all SNPs with info
      N_info = sum((N*INFO)/INFO, na.rm = TRUE),
      N_EAF = sum((N*EAF)/EAF, na.rm = TRUE),
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      B = B / W,
      SE = 1 / sqrt(W),
      dplyr::across(dplyr::any_of(c("INFO")), ~.x / N_info),
      dplyr::across(dplyr::any_of(c("EAF")), ~.x / N_EAF)
    ) |>
    dplyr::select(-W, -N_info, -N_EAF) |>
    dplyr::collect() |>
    dplyr::mutate(P  = stats::pnorm(-abs(B/SE)) *2)

}
