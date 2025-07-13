utils::globalVariables(c("WB2", "n_contributions", "Q", "Q_df"))
#' Improved meta-analysis using tidyGWAS:ed files
#'
#' @param ds a dataset object, see [arrow::open_dataset()]
#' @param min_EAF Filter on minimal EAF (effect allele frequency) value. (0 - 0.5)
#' @param ref Reference genome version to use for reference allele
#' @param chromosomes Which chrosomes to apply meta-analysis across. Default is autosomes: 1:22
#'
#' @returns a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' meta_analyse2(ds, min_EAF = 0.01)
#' }
meta_analyse2 <- function(ds, min_EAF = NULL, ref = c("REF_38", "REF_37"), chromosomes = c(1:22)) {
  ref <- rlang::arg_match(ref)
  schema <- arrow::schema(ds)
  dataset_names <- names(schema)
  chr_type <- schema$CHR
  if(!schema$CHR$type == arrow::int32()) {
    chromosomes <- as.character(chromosomes)
  }

  mandatory <- c(ref, "CHR", "RSID", "EffectAllele", "OtherAllele", "B", "SE")
  all(mandatory %in% dataset_names) ||
    cli::cli_abort(
      "The dataset is missing mandatory columns: {.arg {mandatory[!mandatory %in% dataset_names]}} ",
    )

  # to not break earlier versions of tidyGWAS
  if("indel" %in% dataset_names) {
    ds <- dplyr::filter(ds, is.na(indel) | !indel)
  }

  ds <- dplyr::filter(ds, dplyr::if_all(dplyr::all_of(c("CHR", "RSID", "EffectAllele","OtherAllele")), ~!is.na(.x)))

  if(!is.null(min_EAF)) {
    "EAF" %in% dataset_names ||
    rlang::abort("`min_EAF` is set, but the dataset does not contain an `EAF` column")
    rlang::is_scalar_double(min_EAF) || rlang::abort("`min_EAF` must be a single numeric value.")
    ds <- dplyr::filter(ds, EAF >= min_EAF & EAF <= min_EAF)
  }

  if(!("N" %in% dataset_names)) {
    cli::cli_warn("No N in columns. Cannot average INFO and EAF")
    ds <- dplyr::select(ds, -dplyr::any_of(c("EAF", "INFO")))
  }

  # ----------------------------------------------------------
  # 3. Run the per‑chromosome meta‑analysis in parallel
  # ----------------------------------------------------------
  purrr::map(
    chromosomes,
    \(.x) by_chrom(ds, chrom = .x, ref=ref),
    .progress = list(type = "tasks")
  ) |>
    purrr::list_rbind()
}

by_chrom <- function(ds, chrom, ref) {

  by = c("RSID", "EffectAllele", "OtherAllele")
  cols <- c("B","SE", "N", "CaseN", "ControlN", "EffectiveN","EAF", "INFO", by, "POS_38", "POS_37")

  ds |>
    dplyr::filter(CHR == chrom) |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(c("B", "SE")), ~ is.finite(.x))) |>
    dplyr::rename(REF = !!ref) |>
    align_to_ref() |>
    dplyr::select(dplyr::any_of(cols)) |>
    dplyr::mutate(
      W       = 1 / (SE^2),
      B       = B * W,
      WB2  = (B^2) / W,
      dplyr::across(dplyr::any_of(c("EAF", "INFO")), ~ if_else(is.na(N), 0, N),  .names = "N_{.col}"),
      dplyr::across(dplyr::any_of(c("EAF", "INFO")),  ~ if_else(is.na(N), 0, .x)),
      dplyr::across(dplyr::any_of("EAF"), ~.x * N_EAF),
      dplyr::across(dplyr::any_of("INFO"), ~.x * N_INFO)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      POS_38 = min(POS_38),
      POS_37 = min(POS_37),
      n_contributions = dplyr::n(),
      WB2 = sum(WB2),
      dplyr::across(dplyr::any_of(c("W", "B")), sum),
      dplyr::across(dplyr::any_of(c("EAF", "INFO","N_EAF", "N_INFO","CaseN", "ControlN", "N", "EffectiveN")), ~sum(.x, na.rm=T))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      B   = B / W,
      SE  = 1 / sqrt(W),
      dplyr::across(dplyr::any_of("EAF"), ~.x / N_EAF),
      dplyr::across(dplyr::any_of("INFO"), ~.x / N_INFO),
      CHR = {{chrom}}
    ) |>
    dplyr::select(-dplyr::any_of(c("N_INFO", "N_EAF"))) |>
    dplyr::collect() |>
    dplyr::mutate(
      P       = stats::pnorm(-abs(B / SE)) * 2,
      Q        = WB2 - (B^2) / W,
      Q_df     = pmax(n_contributions - 1L, 0L),
      Q_pval   = stats::pchisq(Q, df = Q_df, lower.tail = FALSE),
      I2       = dplyr::if_else(Q > 0,pmax((Q - Q_df) / Q, 0),0)
    ) |>
    dplyr::select(-"W", "WB2")

}

