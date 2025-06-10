utils::globalVariables(c("MAF"))
#' Repair statistics column in a GWAS summary statistics tibble
#'
#' `repair_stats()` is a collection of functions that can be used to
#' infer missing columns in GWAS summary statistics. The functions are based on
#' functionality found online.
#' @inheritParams tidyGWAS
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' updated <- repair_stats(my_gwas)
#' }
repair_stats <- function(tbl, dbsnp_path, impute_freq = c("None", "EUR", "AMR", "AFR", "SAS", "EAS"), impute_freq_file = NULL, impute_n=FALSE) {

  start_cols <- colnames(tbl)
  impute_freq <- rlang::arg_match(impute_freq)

  # check if allele frequency should be repaired
  if(impute_freq != "None" & is.null(impute_freq_file)) {
    cli::cli_alert_info("Imputing allele frequency using 1000G data, using population: {impute_freq}")
    df_eaf <- arrow::open_dataset(paste0(dbsnp_path, "/EAF_REF_1KG")) |>
      dplyr::filter(.data[["ancestry"]] == impute_freq) |>
      dplyr::select(CHR, POS_38 = POS, EffectAllele, OtherAllele, EAF) |>
      dplyr::collect()
    repair_eaf <- TRUE

  } else if(!is.null(impute_freq_file)) {

    cli::cli_alert_info("Imputing allele frequency using custom frequency file")
    df_eaf <- arrow::read_parquet(impute_freq_file)
    stopifnot(all(c("RSID", "EffectAllele", "OtherAllele", "EAF") %in% colnames(df_eaf)))
    repair_eaf <- TRUE

  } else {

    repair_eaf <- FALSE

  }



  if(repair_eaf & !"EAF" %in% colnames(tbl)) {


    tbl1 <- dplyr::inner_join(tbl, df_eaf, by = c("CHR","POS_38", "EffectAllele", "OtherAllele")) |>
      dplyr::select(c("CHR","POS_38", "EffectAllele", "OtherAllele", "EAF"))

    tbl2 <- dplyr::inner_join(tbl, df_eaf, by = c("CHR","POS_38", "EffectAllele" ="OtherAllele", "OtherAllele" = "EffectAllele")) |>
      dplyr::mutate(EAF = 1-EAF) |>
      dplyr::select(c("CHR","POS_38", "EffectAllele", "OtherAllele", "EAF"))

    final_eaf <- dplyr::bind_rows(tbl1, tbl2)
    cli::cli_alert_success("Was able to add allele frequency for {nrow(final_eaf)} / {nrow(tbl)} SNPs")
    tbl <- dplyr::left_join(tbl, final_eaf, by = c("CHR","POS_38", "EffectAllele", "OtherAllele"))



  }
  if(impute_n & !"N" %in% colnames(tbl) & all(c("EAF", "SE") %in% colnames(tbl))) {
    cli::cli_alert_info("N missing: Calculating N using the formula:  N = 4/( (2 * MAF * (1-MAF)) * SE^2)")
    tbl <- impute_N(tbl)
  }

  if(all(c("B", "SE") %in% colnames(tbl)) & !"Z" %in% colnames(tbl)) {
    cli::cli_alert_info("Z missing: Calculating Z using the formula:  Z = B / SE")
    tbl <- dplyr::mutate(tbl, Z = B / SE)
  }

  if(all(c("B", "P") %in% colnames(tbl)) & !"Z" %in% colnames(tbl)) {
    cli::cli_alert_info("Found B and P but not Z. Imputing Z using:
                        sign(beta) * sqrt(stats::qchisq(pvalue,1,lower=FALSE))")
    tbl <- dplyr::mutate(tbl, Z = z_from_p_b(P, B))
  }

  if(all(c("B", "Z") %in% colnames(tbl)) & !"SE" %in% colnames(tbl)) {
    cli::cli_alert_info("Found B and Z but not SE. Calculating SE using:
                        B / Z = SE")
    tbl <- dplyr::mutate(tbl, SE = B / Z)
  }


  if(all(c("Z", "N", "EAF") %in% colnames(tbl)) & all(!c("B", "SE") %in% colnames(tbl))) {
    cli::cli_alert_info("{.strong Imputing B and SE, using Z, EAF and N }. The B and SE will correspond to a standardized scale, which might not always be the same scale that the original betas was on.")
    # check variable N
    if(length(unique(tbl$N)) == 1) cli::cli_alert_danger("Attempting to repair B and SE from Z,EAF and N. However, N does not not vary across SNPs, indicating it's the study-wide N and not SNP-wise N. This will introduce issues with the conversion")
    tbl <- dplyr::mutate(tbl,
                         B = beta_from_z_eaf_n(Z,EAF,N),
                         SE = se_from_z_eaf_n(Z,EAF,N)
    )
  }
  if("Z" %in% colnames(tbl) & !"P" %in% colnames(tbl)){
    cli::cli_alert("Imputing P based on Z")
    tbl <- dplyr::mutate(tbl,P = stats::pnorm(-abs(tbl$Z)) *2)

  }



  # end message -------------------------------------------------------------

  end_cols <- colnames(tbl)
  new_cols <- end_cols[!end_cols %in% start_cols]
  if(length(new_cols) > 0) {
    cli::cli_h3("Finished repair_stats: ")
    cli::cli_alert_info("Added {length(new_cols)} new columns: {new_cols}")

  }


  # finished - return tibble ---------------------------------------------------
  tbl

}


beta_from_z_eaf_n <- function(Z, EAF, N) {
  # ref https://www.biostars.org/p/319584/
  # and https://www.nature.com/articles/ng.3538
  Z / sqrt((2*EAF*(1 - EAF)) * (N + (Z^2)))
}

se_from_z_eaf_n <- function(Z, EAF, N) {
  # ref https://www.biostars.org/p/319584/
  # https://www.nature.com/articles/ng.3538
  1/sqrt((2*EAF)*(1-(EAF))*(N+(Z^2)))
}

z_from_p_b <- function(pvalue, beta) {
  sign(beta) * sqrt(stats::qchisq(pvalue,1,lower=FALSE))
}

n_from_se_eaf <- function(SE, EAF) {
  # source: https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics
  # convert all estimats to minor allele frequency
  EAF <- dplyr::if_else(EAF > 0.5, 1-EAF, EAF)
  round(4 / ((2 * EAF * (1-EAF)) * (SE^2) ))
}

impute_N <- function(tbl) {
  # source: https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics
  columns <- colnames(tbl)
  stopifnot("EAF" %in% columns)
  stopifnot("SE" %in% columns)

  tbl |>
    dplyr::mutate(
      MAF = dplyr::if_else(EAF > 0.5, 1 - EAF, EAF),
      N = 4/( (2 * MAF * (1-MAF)) * SE^2),
      N = as.integer(N)
    ) |>
    dplyr::select(-MAF)

}
