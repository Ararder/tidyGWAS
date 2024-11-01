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
repair_stats <- function(tbl) {

  start_cols <- colnames(tbl)

  if(all(c("B", "SE") %in% colnames(tbl)) & !"Z" %in% colnames(tbl)) {
    cli::cli_alert_info("Z missing: Calculating Z using the formula:  Z = B / SE")
    tbl <- dplyr::mutate(tbl, Z = B / SE)
  }

  if(all(c("B", "P") %in% colnames(tbl)) & !"Z" %in% colnames(tbl)) {
    cli::cli_alert_info("Found B and P but not Z. Imputing Z using:
                        sign(beta) * sqrt(stats::qchisq(pvalue,1,lower=FALSE))")
    tbl <- dplyr::mutate(tbl, Z = z_from_p_b(P, B))
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
