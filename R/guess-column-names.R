



guess_names <- function(tbl) {

  # find best mapping
  cols <- colnames(tbl)
  scores <- purrr::map_dbl(formats, \(x) count_matches(cols, x))
  best_format_name <- names(scores)[which.max(scores)]
  best_format      <- formats[[best_format_name]]




  success <- best_format[best_format %in% cols]

  cli::cli_h3("The detected format is {.val {best_format_name}}")
  cli::cli_inform("Was able to map {length(success)} out of {length(cols)} columns to tidyGWAS columns")
  for(i in seq_along(success)) {
    cli::cli_alert_success(
      "{names(success)[i]} -> {success[i]}"
    )
  }

  missing <- setdiff(cols, unname(unlist(best_format)))
  # rowid should not be reported as an unmapped column
  missing <- missing[missing != "rowid"]
  if (length(missing)) cli::cli_alert_warning(
    "extra column(s) {.val {missing}} were not mapped to a tidyGWAS column"
  )



  # finished ----------------------------------------------------------------
  dplyr::rename(tbl, !!!success)



}

count_matches <- function(cnames, fmt) {
  fmt_cols <- unname(unlist(fmt))
  sum(fmt_cols %in% cnames)
}






formats <- list(
  SSF  = list(
    CHR = "chromosome",
    POS = "base_pair_location",
    RSID = "rsid",
    EffectAllele = "effect_allele",
    OtherAllele = "other_allele",
    B = "beta",
    SE = "standard_error",
    EAF = "effect_allele_frequency",
    P = "p_value",
    N = "n"
  ),

  REGENIE = list(
    CHR = "CHROM",
    POS = "GENPOS",
    RSID = "ID",
    EffectAllele = "ALLELE1",
    OtherAllele = "ALLELE0",
    EAF = "A1FREQ",
    B = "BETA",
    CaseN = "N_CASES",
    ControlN = "N_CONTROLS",
    INFO = "INFO",
    SE = "SE",
    N = "N"
  ),
  # source: https://yanglab.westlake.edu.cn/software/gcta/index.html#DataResource
  fastGWA = list(
    CHR = "CHR",
    POS = "POS",
    RSID = "SNP",
    EffectAllele = "A1",
    OtherAllele = "A2",
    EAF = "AF1",
    B = "BETA",
    N = "N",
    P = "P",
    SE = "SE"
  ),
  plink = list(
    POS = "BP",
    RSID = "SNP",
    EffectAllele = "A1",
    OtherAllele = "A2",
    EAF = "MAF",
    N = "NMISS"

  ),
  plink_glm = list(
      CHR = "#CHROM",
      RSID = "ID",
      EffectAllele = "A1",
      OtherAllele = "REF",
      EAF = "A1_FREQ",
      SE = "LOG(OR)_SE",
      Z = "Z_STAT",
      N = "OBS_CT"
  )

)




