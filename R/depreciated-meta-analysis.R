utils::globalVariables(
  c("REF", "EA_is_ref", "CHR", "ID", "W", "N_info", "INFO","tmp", "EffectiveN", "N_EAF", "tdg_path", "name")
)


#' Perform meta-analysis of GWAS summary statistics datasets cleaned by tidyGWAS
#' This function is depreciated. Use [meta_analyse()]
#'
#'
#' @param dset an [arrow::open_dataset()] object
#' @param by a character vector of column names to group by. Default is c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele")
#' @param ref either "REF_37" or "REF_38", depending on which column you want to use to standardize reference allele
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' dset <- arrow::open_dataset("path_to/sumstats/")
#' res <- meta_analyze(dset)
#' }
#'
deprec_meta_analyze <- function(dset, by = c("CHR", "POS_37", "RSID", "EffectAllele", "OtherAllele"), ref = c("REF_37", "REF_38")) {
  ref <- rlang::arg_match(ref)

  purrr::map(c(1:22), \(chrom) meta_analyze_by_chrom(dset, chrom = chrom, by = by, ref = ref)) |>
    purrr::list_rbind()

}



meta_analyze_by_chrom <- function(dset, chrom, by, ref) {
  stats <- c("B", "SE", "EAF", "N", "CaseN", "ControlN","INFO", "EffectiveN")
  cols <- c(by, stats)

  q1 <- dplyr::filter(dset, CHR == {{ chrom }}) |>
   dplyr::rename(REF = dplyr::all_of(ref))

  if("indel" %in% names(arrow::schema(dset))) {
     q1 <- dplyr::filter(q1, is.na(indel) | !indel)
  }

  dset |>
    dplyr::filter(CHR == {{ chrom }}) |>
    dplyr::rename(REF = !!ref) |>
    align_to_ref() |>
    dplyr::select(dplyr::any_of(cols)) |>
    dplyr::mutate(
      W = 1 / (SE^2),
      B = B*W,
      dplyr::across(dplyr::any_of(c("EAF", "INFO")), ~.x * N)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    handle_info_eaf() |>
    dplyr::select(-dplyr::any_of(c("W", "N_info", "N_EAF"))) |>
    dplyr::collect() |>
    dplyr::mutate(P  = stats::pnorm(-abs(B/SE)) *2)

}




handle_info_eaf <- function(query) {
  if(all(c("INFO", "EAF") %in% names(arrow::schema(query)))) {
    query |>
      dplyr::summarise(
        N_info = sum( N * (INFO/INFO), na.rm = TRUE),
        N_EAF = sum(  N * (EAF/EAF), na.rm = TRUE),
        n_contributions = dplyr::n(),
        dplyr::across(dplyr::any_of(c("W", "B")), sum),
        dplyr::across(dplyr::all_of(c("EAF", "INFO")), ~sum(.x, na.rm=T)),
        dplyr::across(dplyr::any_of(c("CaseN", "ControlN", "N", "EffectiveN")), ~sum(.x, na.rm=T)),
        # calculate the total sample size for all SNPs with info
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        B = B / W,
        SE = 1 / sqrt(W),
        INFO = INFO / N_info,
        EAF = EAF / N_EAF
      )


  } else if("INFO" %in% names(arrow::schema(query))) {
    query |>
      dplyr::summarise(
        N_info = sum( N * (INFO/INFO), na.rm = TRUE),
        n_contributions = dplyr::n(),
        dplyr::across(dplyr::any_of(c("W", "B")), sum),
        dplyr::across(dplyr::all_of(c("INFO")), ~sum(.x, na.rm=T)),
        dplyr::across(dplyr::any_of(c("CaseN", "ControlN", "N", "EffectiveN")), ~sum(.x, na.rm=T)),
        # calculate the total sample size for all SNPs with info
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        B = B / W,
        SE = 1 / sqrt(W),
        INFO = INFO / N_info,
      )



  } else if("EAF" %in% names(arrow::schema(query))) {
    query |>
      dplyr::summarise(
        N_EAF = sum( N*(EAF/EAF), na.rm = TRUE),
        n_contributions = dplyr::n(),
        dplyr::across(dplyr::any_of(c("W", "B")), sum),
        dplyr::across(dplyr::all_of(c("EAF")), ~sum(.x, na.rm=T)),
        dplyr::across(dplyr::any_of(c("CaseN", "ControlN", "N", "EffectiveN")), ~sum(.x, na.rm=T)),
        # calculate the total sample size for all SNPs with info
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        B = B / W,
        SE = 1 / sqrt(W),
        EAF = EAF / N_EAF
      )


  } else {
    query |>
      dplyr::summarise(
        n_contributions = dplyr::n(),
        dplyr::across(dplyr::any_of(c("W", "B")), sum),
        dplyr::across(dplyr::any_of(c("CaseN", "ControlN", "N", "EffectiveN")), ~sum(.x, na.rm=T))
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        B = B / W,
        SE = 1 / sqrt(W)
      )


  }

}


