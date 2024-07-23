utils::globalVariables(c(
  "RSID", "POS", "CHR", "EffectAllele", "OtherAllele",
  "seqnames", "pos", "RefSNP_id", "ref_allele", "alt_alleles","REF_38",
  "a1_is_ref", "a2_is_ref","incompat_alleles", "rowid", "RSID.x", "POS_38", ":=", "RSID.y","convert_p"
  ))



#' Infer what genome build a GWAS summary statistics file is on.
#' @inheritParams tidyGWAS
#' @param n_snps number of snps to check CHR and POS for
#'
#' @return either "37" or "38"
#' @export
#'
#' @examples \dontrun{
#' genome_build <- infer_build(gwas_sumstats)
#' }
infer_build <- function(tbl, dbsnp_path, n_snps = 10000) {
  stopifnot("Need 'CHR' and 'POS' in tbl" = all(c("CHR", "POS") %in% colnames(tbl)))
  cli::cli_alert_info("Inferring build by matching {n_snps} rows to GRCh37 and GRCh38")

  subset <- dplyr::slice_sample(tbl, n = {{ n_snps }})
  b38 <- map_to_dbsnp(tbl = subset, build = "38", dbsnp_path = dbsnp_path)
  b37 <- map_to_dbsnp(tbl = subset, build = "37", dbsnp_path = dbsnp_path)

  if(nrow(b37) > nrow(b38)) build <- "37" else build <- "38"
  cli::cli_inform("{nrow(b38)} snps matched GRCh38, {nrow(b37)} for GRCh37, inferring build to be {build}")

  build


}


# workhorse of performing joins between the ~1 billion rows of dbSNP and summary statistics to be clenaed
map_to_dbsnp <- function(tbl, dbsnp_path, build = c("NA", "37", "38")) {


  rlang::check_required(tbl)
  rlang::check_required(dbsnp_path)
  build <- rlang::arg_match(build)


  if(!dplyr::is.tbl(tbl)) tbl <- dplyr::tibble(tbl)
  stopifnot("Cannot map 0 length tibble to dbsnp" = nrow(tbl) > 0)




  dset <- arrow::open_dataset(paste0(dbsnp_path, "/v155"))


  # join with dbsnp ----------------------------------------------------------
  if(build == "NA") {
    stopifnot("When mapping with RSID, RSID needs to present" = "RSID" %in% colnames(tbl))
    stopifnot(is.character(tbl$RSID))

    # RSID is stored as integer with "rs" suffix removed for better performance
    tbl$RSID <- as.integer(stringr::str_sub(tbl$RSID, start = 3))

    # batching by CHR gives SIGNIFICANTLY better performance even though we check ALL rows against each chromosome
    chrom <- c(1:22, "X", "Y", "MT")

    res <-
      purrr::map(chrom, \(chrom) dplyr::filter(dset, CHR == {{ chrom }} & RSID %in% tbl$RSID) |> dplyr::collect()) |>
      purrr::list_rbind()


  } else {

    stopifnot("'CHR' and 'POS' need to be present in tbl" =  all(c("CHR", "POS") %in% colnames(tbl)))
    tbl$CHR <- as.character(tbl$CHR)
    tbl$POS <- as.integer(tbl$POS)

    # we remove variants that do not have the same chromosome across builds - therefore we can use CHR_38 regardless of build
    if(build == "37")
      res <-
      split(tbl, tbl$CHR) |>
      purrr::imap(\(df, chrom) dplyr::filter(dset, CHR == {{ chrom }} & POS_37 %in% df$POS) |> dplyr::collect()) |>
      purrr::list_rbind()
    else {
      res <-
        split(tbl, tbl$CHR) |>
        purrr::imap(\(df, chrom) dplyr::filter(dset, CHR == {{ chrom }} & POS_38 %in% df$POS) |> dplyr::collect()) |>
        purrr::list_rbind()

    }

  }



  res |>
    dplyr::mutate(RSID = stringr::str_c("rs", RSID)) |>
    # ensure correct types
    dplyr::mutate(
      dplyr::across(c("CHR", "REF_38", "ALT_38", "REF_37", "ALT_37"), as.character),
      POS_37 = as.integer(POS_37),
      POS_38 = as.integer(POS_38)
    )



}








repair_ids <- function(tbl, repair = c("rsid", "pos"), build = c("NA", "37", "38"), dbsnp_path) {
  repair <- rlang::arg_match(repair)
  build <- rlang::arg_match(build)
  rlang::check_required(dbsnp_path)

  if(repair == "rsid") {
    tbl <- dplyr::select(tbl, -dplyr::any_of("RSID"))
    if(build == "NA") build <- infer_build(tbl, dbsnp_path = dbsnp_path)
    merge_cols <- c("CHR", paste0("POS_", build))
    cols_for_allele_check <- c("CHR", paste0("POS_", build), "RSID", "REF_37", "REF_38")
    dbsnp <- map_to_dbsnp(tbl, build = build, dbsnp_path = dbsnp_path)

    new_name <- paste0("POS_", build)
    tbl <- dplyr::rename(tbl, {{ new_name }} := "POS")
  } else {
    merge_cols <- c("RSID")
    # if matching by RSID, there is no concept of genome build. However,
    # there's a subset of entries where dbSNP ref allele differs by build.
    # most GWASes are still on GRCh37, so we default to that
    build <- "37"
    tbl <- dplyr::select(tbl, -dplyr::any_of(c("CHR", "POS")))
    if(!"rowid" %in% colnames(tbl)) tbl$rowid <- 1:nrow(tbl)
    dbsnp <- map_to_dbsnp(tbl = tbl, dbsnp_path = dbsnp_path)

  }

  # append CHR/POS or RSID
  tbl <-
    dplyr::left_join(
      tbl,
      dplyr::select(dbsnp, dplyr::all_of(c("CHR", "POS_38", "POS_37", "RSID", "REF_37", "REF_38"))),
      by = merge_cols
    ) |>
    # merge_cols can be either length 1 or 2 depending on matching by rsid or CHR-POS
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(REF_38), TRUE, FALSE))




  flags <- check_incompat_alleles(tbl, dbsnp_df = dbsnp, build = build)


  # merge in flag
  dplyr::left_join(tbl, flags, by = "rowid")

}


check_incompat_alleles <- function(tbl, dbsnp_df, build = c("37", "38")) {

  build <- rlang::arg_match(build)
  req_cols1 <- c("rowid", "EffectAllele", "OtherAllele", "RSID")
  stopifnot("Missing mandatory columns from tbl" = all(req_cols1 %in% colnames(tbl)))
  req_cols <- c("RSID", "CHR","POS_37","POS_38", "REF_38", "REF_37", "ALT_37", "ALT_38")
  stopifnot("Missing mandatory columns from dbsnp_df" = all(req_cols %in% colnames(dbsnp_df)))

  if(build == "37") {
    dbsnp_cols <- c("CHR", "POS" = "POS_37", "RSID", "ref_allele" = "REF_37", "alt_alleles" = "ALT_37")
  } else {
    dbsnp_cols <- c("CHR", "POS" = "POS_38", "RSID", "ref_allele" = "REF_38", "alt_alleles" = "ALT_38")
  }


  # correct cols
  tbl <- dplyr::select(tbl, dplyr::all_of(c("rowid", "RSID", "EffectAllele", "OtherAllele")))

  # handle multi-allelics
  dbsnp_df <- dplyr::select(dbsnp_df, dplyr::all_of(dbsnp_cols)) |>
    # multi-allelic SNPs into separate rows
    tidyr::separate_longer_delim(alt_alleles, delim =",")



  # check that EffectAllele and OtherAllele is either REF/ALT, or ALT/REF
  no_strand_flip <- flag_incompat_alleles(tbl, dbsnp_df)
  failed <- dplyr::anti_join(tbl, no_strand_flip, by = "rowid")
  strand_flip <- flag_incompat_alleles(strand_flip(dplyr::anti_join(tbl, no_strand_flip, by = "rowid")), dbsnp_df)
  merged <- dplyr::bind_rows(no_strand_flip, strand_flip)

  # finished ----------------------------------------------------------------
  dplyr::mutate(tbl, incompat_alleles = dplyr::if_else(rowid %in% merged$rowid, FALSE, TRUE)) |>
    dplyr::select(dplyr::all_of(c("rowid", "incompat_alleles")))


}

# -------------------------------------------------------------------------


flag_incompat_alleles <- function(tbl, dbsnp_df) {
  # if input contains
  dplyr::inner_join(tbl, dbsnp_df, by = "RSID", relationship = "many-to-many") |>
    dplyr::mutate(
      a1_is_ref = dplyr::if_else(EffectAllele == ref_allele & OtherAllele == alt_alleles, TRUE, FALSE),
      a2_is_ref = dplyr::if_else(OtherAllele == ref_allele & EffectAllele == alt_alleles, TRUE, FALSE),
    ) |>
    dplyr::mutate(incompat_alleles = dplyr::if_else(a1_is_ref | a2_is_ref, FALSE, TRUE)) |>
    dplyr::select(-a1_is_ref, -a2_is_ref, -ref_allele, -alt_alleles) |>
    dplyr::filter(!incompat_alleles)
}






