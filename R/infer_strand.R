#' utils::globalVariables(c("EffectAllele", "OtherAllele", "EAF", "B", "EAF_ref", "diff_direct", "diff_flipped"))
#'
#' #' Infer and correct strand for palindromic SNPs
#' #'
#' #' This function first aligns the summary statistics to a reference panel using
#' #' `align_to_ref`. Then, for the remaining ambiguous palindromic SNPs, it
#' #' compares the allele frequencies against a reference panel to infer and
#' #' correct likely strand flips.
#' #'
#' #' @param tbl A data.frame or tibble with GWAS summary statistics. Must contain
#' #'   `rowid`, `CHR`, `POS`, `EffectAllele`, `OtherAllele`, `EAF`, and `B`.
#' #' @param ref_freq A data.frame or tibble with reference allele frequencies.
#' #'   Must contain `CHR`, `POS`, `EffectAllele`, `OtherAllele`, and `EAF`.
#' #' @param ref_build The reference build to use for the initial alignment,
#' #'   either "REF_37" or "REF_38".
#' #'
#' #' @return A tibble with corrected strands, effects, and a new `inferred_flip`
#' #'   column indicating which rows were corrected.
#' #' @export
#' infer_palindromic_strand <- function(tbl, ref_freq, ref_build = c("REF_38", "REF_37")) {
#'   ref_build <- rlang::arg_match(ref_build)
#'
#'   # 1. Standardize alleles against the reference panel first
#'   cli::cli_alert_info("Step 1: Aligning all SNPs to the reference panel.")
#'   aligned_tbl <- align_to_ref(tbl, ref = ref_build)
#'
#'   # 2. Isolate palindromic SNPs for ambiguity resolution
#'   palindromic_snps <- aligned_tbl |>
#'     dplyr::filter(
#'       (EffectAllele == "A" & OtherAllele == "T") |
#'       (EffectAllele == "T" & OtherAllele == "A") |
#'       (EffectAllele == "C" & OtherAllele == "G") |
#'       (EffectAllele == "G" & OtherAllele == "C")
#'     )
#'
#'   if (nrow(palindromic_snps) == 0) {
#'     cli::cli_alert_info("No palindromic SNPs found after alignment. Returning.")
#'     return(dplyr::mutate(aligned_tbl, inferred_flip = FALSE))
#'   }
#'
#'   cli::cli_alert_info("Step 2: Found {nrow(palindromic_snps)} palindromic SNPs to check.")
#'
#'   # 3. Merge and compare frequency scenarios
#'   merged_data <- dplyr::inner_join(
#'     palindromic_snps,
#'     ref_freq,
#'     by = c("CHR", "POS", "EffectAllele", "OtherAllele"),
#'     suffix = c("", "_ref")
#'   )
#'
#'   if (nrow(merged_data) == 0) {
#'     cli::cli_alert_info("No matching palindromic SNPs in the reference. No flips inferred.")
#'     return(dplyr::mutate(aligned_tbl, inferred_flip = FALSE))
#'   }
#'
#'   # 4. Infer flips based on which scenario is more likely
#'   inferred_snps <- merged_data |>
#'     dplyr::mutate(
#'       diff_direct = abs(EAF - EAF_ref),
#'       diff_flipped = abs(EAF - (1 - EAF_ref)),
#'       # Decision: A flip is inferred if the flipped frequency is a much better match
#'       inferred_flip = diff_flipped < 0.15 & diff_direct > 0.85
#'     )
#'
#'   snps_to_correct <- dplyr::filter(inferred_snps, inferred_flip)
#'
#'   if (nrow(snps_to_correct) > 0) {
#'     cli::cli_alert_info("Step 3: Inferred and corrected {nrow(snps_to_correct)} strand flips.")
#'
#'     # 5. Correct the SNPs that were flagged
#'     corrected_rows <- snps_to_correct |>
#'       dplyr::mutate(
#'         EffectAllele_new = OtherAllele,
#'         OtherAllele_new = EffectAllele,
#'         EffectAllele = EffectAllele_new,
#'         OtherAllele = OtherAllele_new,
#'         B = -B,
#'         EAF = 1 - EAF
#'       ) |>
#'       dplyr::select(dplyr::all_of(colnames(tbl)), inferred_flip)
#'
#'     # Update the main table and add the tracking column
#'     final_tbl <- dplyr::rows_update(aligned_tbl, corrected_rows, by = "rowid") |>
#'       dplyr::left_join(dplyr::select(corrected_rows, rowid, inferred_flip), by = "rowid") |>
#'       dplyr::mutate(inferred_flip = tidyr::replace_na(inferred_flip, FALSE))
#'
#'   } else {
#'     cli::cli_alert_info("Step 3: No strand flips inferred based on frequency comparison.")
#'     final_tbl <- dplyr::mutate(aligned_tbl, inferred_flip = FALSE)
#'   }
#'
#'   return(final_tbl)
#' }
