validate_sumstats_only_rsid <- function(with_chr_pos=TRUE) {

  if(!with_chr_pos) {

  cli::cli_alert_info(
    "These rows were removed when validating columns.
    `tidyGWAS()` detected only an RSID column and no CHR POS columns to create a variant id.
    "
    )
  } else {
    cli::cli_alert_info(
      "These rows were removed when validating columns.
    `tidyGWAS()` detected only an RSID column and no CHR POS columns to create a variant id.
    During this process, a subset of rows were detected as having CHR:POS or CHR:POS:REF:ALT
    in the RSID column. From this subset of rows the following rows were removed.
    "
      )
  }
}
