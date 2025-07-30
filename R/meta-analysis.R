utils::globalVariables(c("WB2", "n_contributions", "Q", "Q_df", "B_w", "Q_pval"))

#' Improved meta-analysis using tidyGWAS:ed files
#'
#'  @description
#'  [meta_analyse()] will:
#' - flip the effect allele to be the reference allele (using either GRCh37 or GRCh38), control with `ref`
#' - variant id is constructed using RSID:EffectAllele:OtherAllele (Which is now RSID:REF:ALT)
#' - EAF (allele frequency) and INFO (imputation quality) are weighted by sample size, if present
#' - CaseN, N, ControlN and EffectiveN are all summed and carried forward
#'
#'
#' @param ds a dataset object, see [arrow::open_dataset()]
#' @param min_EAF Filter on minimal EAF (effect allele frequency) value. (0 - 0.5)
#' @param ref Reference genome version to use for reference allele
#' @param chromosomes Which chrosomes to apply meta-analysis across. Default is autosomes: 1:22
#' @param safe_mode Apply additional filters to ensure no missing/inf/NaN across key columns?
#'  This is set to FALSE by default because these filters should already be
#'  guaranteed with the use of [tidyGWAS()]
#'
#' @returns a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' meta_analyse2(ds, min_EAF = 0.01)
#' }
meta_analyse <- function(ds, chromosomes = c(1:22), min_EAF = NULL, ref = c("REF_38", "REF_37"), safe_mode = FALSE) {
  ref <- rlang::arg_match(ref)
  rlang::is_scalar_logical(safe_mode) || cli::cli_abort("`safe_mode` must be a single logical value.")
  schema <- arrow::schema(ds)
  dataset_names <- names(schema)

  mandatory <- c(ref, "CHR", "RSID","POS_37","POS_38", "EffectAllele", "OtherAllele", "B", "SE")
  all(mandatory %in% dataset_names) ||
    cli::cli_abort(
      "The dataset is missing mandatory columns: {.arg {mandatory[!mandatory %in% dataset_names]}} ",
    )

  chr_type <- schema$CHR
  if(!schema$CHR$type == arrow::int32()) {
    chromosomes <- as.character(chromosomes)
  }

  # # variant identity columns cannot be NA
  if(safe_mode) {
    ds <- dplyr::filter(ds, dplyr::if_all(dplyr::all_of(c("CHR", "RSID", "EffectAllele","OtherAllele")), ~!is.na(.x))) |>
      dplyr::filter(dplyr::if_all(dplyr::all_of(c("B", "SE")), ~ is.finite(.x)))
  }


  if(!is.null(min_EAF)) {
    "EAF" %in% dataset_names ||
    rlang::abort("`min_EAF` is set, but the dataset does not contain an `EAF` column")
    rlang::is_scalar_double(min_EAF) || rlang::abort("`min_EAF` must be a single numeric value.")
    ds <- dplyr::filter(ds, EAF >= min_EAF & EAF <= (1-min_EAF))
  }

  if(!("N" %in% dataset_names)) {
    cli::cli_alert_warning("No N in columns. Cannot average INFO and EAF")
    ds <- dplyr::select(ds, -dplyr::any_of(c("EAF", "INFO")))
  }


  # loop --------------------------------------------------------------------

  purrr::map(
    chromosomes,
    \(.x) by_chrom(ds, chrom = .x, ref = ref),
    .progress = list(type = "tasks", name = "Meta-analyzing per chromosome...")
  ) |>
    purrr::list_rbind()


  # -------------------------------------------------------------------------


}




by_chrom <- function(ds, chrom, ref) {

  by = c("RSID", "EffectAllele", "OtherAllele")
  cols <- c("B","SE", "N", "CaseN", "ControlN", "EffectiveN","EAF", "INFO", by, "POS_38", "POS_37")

  ds |>
    dplyr::filter(CHR == chrom) |>
    align_to_ref(dset = _, ref = ref) |>
    dplyr::select(dplyr::any_of(cols)) |>

    # here is the core meta-analysis logic
    # -------------------------------------------------------------------------
    dplyr::mutate(
      W       = 1 / (SE^2),
      B       = B * W,
      WB2  = (B^2) / W,
      dplyr::across(dplyr::any_of(c("EAF", "INFO")), ~.x * N)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_contributions = dplyr::n(),
      # denominator for INFO and N has to handle potential missing values in N.
      # if EAF is missing but N is present, cannot include N in denominator
      dplyr::across(dplyr::any_of(c("EAF")), ~  sum(dplyr::if_else(is.na(.x), 0, N), na.rm = TRUE),  .names = "N_{.col}"),
      dplyr::across(dplyr::any_of(c("INFO")), ~ sum(dplyr::if_else(is.na(.x), 0, N), na.rm = TRUE),  .names = "N_{.col}"),
      dplyr::across(dplyr::any_of(c("POS_38", "POS_37")), ~ min(.x)),
      # yes, this needs to be outside of the across
      WB2 = sum(WB2, na.rm = TRUE),
      dplyr::across(dplyr::any_of(c("W", "B")), sum),
      dplyr::across(dplyr::any_of(c("EAF", "INFO","CaseN", "ControlN", "N", "EffectiveN")), ~ sum(.x, na.rm=TRUE))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      B_w = B,
      B   = B / W,
      SE  = 1 / sqrt(W),
      dplyr::across(dplyr::any_of("EAF"), ~.x / N_EAF),
      dplyr::across(dplyr::any_of("INFO"), ~.x / N_INFO),
      CHR = {{ chrom }}
    ) |>

    # -------------------------------------------------------------------------

    dplyr::select(-dplyr::any_of(c("N_INFO", "N_EAF"))) |>
    dplyr::collect() |>
    dplyr::mutate(
      P       = stats::pnorm(-abs(B / SE)) * 2,
      Q        = WB2 - (B_w^2) / W,
      Q_df     = pmax(n_contributions - 1L, 0L),
      Q_pval   = stats::pchisq(Q, df = Q_df, lower.tail = FALSE),
      Q_pval   = dplyr::if_else(Q_df == 0, 1, Q_pval),
      I2       = dplyr::if_else(Q > 0,pmax((Q - Q_df) / Q, 0),0)
    ) |>
    dplyr::select(-c("W","WB2", "B_w"))

}

#' Align EffectAllele to always be the reference genome allele
#' @description
#' [align_to_ref()] will:
#' 1. Filter any variants where EffectAllele or OtherAllele is not the reference genome allele
#' 2. Flip EffectAllele such that it is always the reference genome allele
#' 3. If OtherAllele is the reference genome allele, direction of B,Z and EAF are flipped.
#'
#'
#' @param dset object created by [arrow::open_dataset()] or [dplyr::tibble()]
#' @param ref Reference genome allele to align to. Varies in ~0.5% of locations
#'
#' @returns a [dplyr::tibble()] or arrow query depending on whether dset is a tibble or arrow dataset
#' @export
#'
#' @examples \dontrun{
#' align_to_ref(tidygwas_df, ref = "REF_38")
#' }
align_to_ref <- function(dset, ref = c("REF_38", "REF_37")) {
  ref <- rlang::arg_match(ref)
  # EffectAllele is harmonized to always be the reference allele
  dset |>
    # This can happen if meta-analyzing datasets that were originally on different genome builds
    dplyr::filter(EffectAllele == .data[[ref]] | OtherAllele == .data[[ref]]) |>
    dplyr::mutate(
      EA_is_ref = dplyr::if_else(EffectAllele == .data[[ref]], TRUE,FALSE),
      tmp = EffectAllele,
      EffectAllele = dplyr::if_else(EA_is_ref, EffectAllele, OtherAllele),
      OtherAllele = dplyr::if_else(EA_is_ref, OtherAllele, tmp),
      B = dplyr::if_else(EA_is_ref, B, B*-1),
      dplyr::across(dplyr::any_of("Z"),  ~ dplyr::if_else(EA_is_ref, .x, .x * -1)),
      dplyr::across(dplyr::any_of(c("EAF")), ~dplyr::if_else(EA_is_ref, .x, 1-.x))
    ) |>
    dplyr::select(-dplyr::all_of(c("EA_is_ref", "tmp")))

}



#' Create a data lake in hivestyle format
#'
#' @param dir a directory containing tidyGWAS cleaned summary statistics.
#'  Each directory should contain a tidyGWAS_hivestyle directory.
#' @param lake a directory to create the data lake in.
#'
#' @returns NULL
#' @export
#'
#' @examples \dontrun{
#' # create_lake("/path/to/dir", "/path/to/lake")
#' }
create_lake <- function(dir, lake) {
  paths <- fs::dir_ls(dir, type = "dir", glob = "*tidyGWAS_hivestyle", recurse = 1)
  fs::dir_create(lake)

  tbl <- dplyr::tibble(
    tdg_path = paths,
    name = fs::path_dir(tdg_path) |> fs::path_file(),
    new_path = fs::path(lake, paste0("dataset_name=",name))
  )

  fs::link_create(tbl$tdg_path, tbl$new_path)
}

