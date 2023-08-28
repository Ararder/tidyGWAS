utils::globalVariables(c(
  "RSID", "POS", "CHR", "EffectAllele", "OtherAllele",
  "seqnames", "pos", "RefSNP_id", "ref_allele", "alt_alleles",
  "a1_is_ref", "a2_is_ref","incompat_alleles", "rowid", "RSID.x"
  ))

# https://www.ncbi.nlm.nih.gov/snp/docs/rs_multi_mapping/

# -------------------------------------------------------------------------


#' Compare CHR, POS and RSID with dbSNP reference data
#'
#' @param tbl a [dplyr::tibble()], formated with [tidyGWAS_columns()]
#' @inheritParams tidyGWAS
#' @return a [dplyr::tibble()] with columns `CHR`, `POS`, `POS_37` and `CHR_37` added
#' @export
#'
#' @examples \dontrun{
#' gwas <- tidyGWAS::test_file
#' # make_callback can be passed a filepath to write out removed rows to.
#' callback <- make_callback("~/output_folder/verify_chr_pos_rsid_removed_rows.tsv")
#' verify_chr_pos_rsid(gwas, bs, build = 37)
#' }
verify_chr_pos_rsid <- function(tbl, build = c("NA", "37", "38"), dbsnp_path, add_missing_build =TRUE) {

  # start -------------------------------------------------------------------
  if(!"rowid" %in% colnames(tbl)) tbl$rowid <- 1:nrow(tbl)
  build = rlang::arg_match(build)
  if(build == "NA") build <- infer_build(tbl, dbsnp_path = dbsnp_path)


  # -------------------------------------------------------------------------
  # merge first with CHR:POS, then attempt to merge with RSID
  # if any rows fail CHR:POS merge.



  # get dbsnp data using CHR and POS
  dbsnp_ref <- map_to_dbsnp(tbl, build = build, by = "chr:pos", dbsnp_path = dbsnp_path)

  # check for rows that doesn't match
  chr_pos_not_in_dbsnp <- dplyr::anti_join(tbl, dbsnp_ref, by = c("CHR", "POS"))

  # update with RSID from dbsnp
  updated <- dplyr::inner_join(dplyr::select(tbl, -RSID), dplyr::select(dbsnp_ref, -alt_alleles), by = c("CHR", "POS"))

  # check how many rows had RSIDs different from what was found in dbSNP
  rsid_was_updated <- dplyr::inner_join(tbl, dbsnp_ref, by = c("CHR", "POS")) |>
    dplyr::filter(RSID.x != RSID.y)

  # if any rows could not be found with CHR:POS mapping, try mapping with RSID
  if(nrow(chr_pos_not_in_dbsnp > 0)) {

    dbsnp_ref_2 <- map_to_dbsnp(chr_pos_not_in_dbsnp, build = build, by = "rsid", dbsnp_path = dbsnp_path)

    updated_2 <- dplyr::inner_join(
      dplyr::select(chr_pos_not_in_dbsnp, -CHR, -POS),
      dplyr::select(dbsnp_ref_2, -alt_alleles),
      by = "RSID"
    )
  } else {
    updated_2 <- dplyr::tibble(RSID ="tmp", .rows = 0)
    dbsnp_ref_2 <- dplyr::tibble(RSID ="tmp", .rows = 0)
  }
  # merge the two dbsnp tables
  dbsnp_ref <- dplyr::bind_rows(dbsnp_ref, dbsnp_ref_2) |>
    dplyr::distinct()

  # merge the two rounds of mapping
  updated <- dplyr::bind_rows(updated, updated_2)

  rows_with_changed_values <- dplyr::bind_rows(rsid_was_updated, updated_2)
  removed <- dplyr::filter(tbl, !rowid %in% updated$rowid)



  # -------------------------------------------------------------------------




  if(nrow(rows_with_changed_values) > 0) {
    cli::cli_alert_info("A total of {nrow(rows_with_changed_values)} rows had incompatible CHR:POS and RSID data. These have been updated")
  }


  # add flag that EffectAllele/OtherAllele is incompatible with REF/ALT in dbSNP
  flags <- check_incompat_alleles(updated, dbsnp_df = dbsnp_ref) |>
    dplyr::select(dplyr::all_of(c("rowid", "incompat_alleles")))



  tbl <-
    dplyr::select(tbl, rowid) |>
    dplyr::left_join(updated, by = "rowid") |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(CHR), TRUE, FALSE)) |>
    dplyr::left_join(flags, by = "rowid")


  #3) add CHR and POS from remaining build ------------------------------------
  if(isTRUE(add_missing_build)) {
    tbl <- add_missing_build(tbl, ifelse(build == "37", "38", "37"), dbsnp_path = dbsnp_path)
  }

  tbl
}


#' Use CHR and POS to get RSID from dbSNP 155
#'
#' This function assumes tidyGWAS column names, see [tidyGWAS_columns()]
#'
#' @inheritParams verify_chr_pos_rsid
#'
#' @return a [dplyr::tibble()] with columns `CHR`, `POS`, `POS_37` and `CHR_37` added
#' @export
#'
#' @examples \dontrun{
#' sumstat_df <- repair_rsid(sumstat, bsgenome_list)
#' }
#'
repair_rsid <- function(tbl, build = c("NA", "37", "38"), dbsnp_path, add_missing_build=TRUE){
  build = rlang::arg_match(build)
  tbl <- dplyr::select(tbl, -dplyr::any_of("RSID"))
  if(!"rowid" %in% colnames(tbl)) tbl$rowid <- 1:nrow(tbl)
  if(build == "NA") build <- infer_build(tbl, dbsnp_path = dbsnp_path)


  # Get RSID using CHR:POS
  dbsnp_ref <- map_to_dbsnp(tbl, build = build, by = "chr:pos", dbsnp_path = dbsnp_path)

  # Add RSID, and flag any rows without a matching dbSNP entry
  tbl <-
    dplyr::left_join(tbl, dplyr::select(dbsnp_ref, CHR,POS, RSID, ref_allele), by = c("CHR","POS")) |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(RSID), TRUE, FALSE))

  # identify rows with incompatible alleles
  flags <- dplyr::select(check_incompat_alleles(tbl, dbsnp_ref), dplyr::all_of(c("rowid", "incompat_alleles")))

  # merge in flag
  tbl <- dplyr::left_join(tbl, flags, by = "rowid")

  # add missing build
  if(isTRUE(add_missing_build)) {
    tbl <- add_missing_build(tbl, ifelse(build == "38", "37", "38"), dbsnp_path = dbsnp_path)
  }

  tbl

}





#' Get CHR and POS using RSID
#'
#' @inheritParams verify_chr_pos_rsid
#' @return a [dplyr::tibble()] with columns `CHR`, `POS`, `POS_37` and `CHR_37` added
#' @export
#'
#' @examples \dontrun{
#' sumstat_df <- repair_chr_pos(sumstat)
#' }
#'
repair_chr_pos <- function(tbl, dbsnp_path, add_missing_build=TRUE) {
  tbl <- dplyr::select(tbl, -dplyr::any_of(c("CHR", "POS")))
  if(!"rowid" %in% colnames(tbl)) tbl$rowid <- 1:nrow(tbl)


  # start -------------------------------------------------------------------------

  # use RSID to get CHR and POS on GRCh38
  dbsnp_38 <- map_to_dbsnp(tbl = tbl, build = "38", by = "rsid", dbsnp_path = dbsnp_path)

  # add CHR and POS, and flag that indicates if a match was found in dbsnp
  tbl <- dplyr::left_join(tbl, dplyr::select(dbsnp_38, CHR, POS, RSID, ref_allele), by = "RSID") |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(CHR), TRUE, FALSE))

  # add flags to see if REF/ALT is consistent with EA/OA
  flags <- check_incompat_alleles(tbl, dbsnp_df = dbsnp_38) |>
    dplyr::select(rowid, incompat_alleles)


  # add GRCh37 CHR POS columns
  if(isTRUE(add_missing_build)) {
    tbl <- add_missing_build(tbl, missing_build = "37", dbsnp_path = dbsnp_path)
  }
  # return ----------------------------------------------------------------

  dplyr::left_join(tbl, flags, by = "rowid")


}


# -------------------------------------------------------------------------


add_missing_build <- function(sumstat, missing_build = c("37", "38"), dbsnp_path) {
  missing_build = rlang::arg_match(missing_build)

  dbsnp <- map_to_dbsnp(sumstat, build = missing_build, by = "rsid", dbsnp_path = dbsnp_path) |>
    dplyr::select(-dplyr::all_of(c("ref_allele", "alt_alleles"))) |>
    # remove multi-allelic SNPs
    dplyr::distinct(CHR, POS, RSID)


  if(missing_build == 38) {
    sumstat <- dplyr::rename(sumstat, CHR_37 = CHR, POS_37 = POS)
  } else {
    dbsnp <- dplyr::rename(dbsnp, CHR_37 = CHR, POS_37 = POS)
  }

  dplyr::left_join(sumstat, dbsnp, by = "RSID")

}

# -------------------------------------------------------------------------



#' Infer what genome build a GWAS summary statistics file is on.
#' @inheritParams verify_chr_pos_rsid
#' @param n_snps number of snps to check CHR and POS for
#'
#' @return either "37" or "38"
#' @export
#'
#' @examples \dontrun{
#' genome_build <- infer_build(gwas_sumstats)
#' }
infer_build <- function(tbl, n_snps = 10000, dbsnp_path) {
  cli::cli_alert_info("Inferring build by matching {n_snps} rows to GRCh37 and GRCh38")
  stopifnot("Need 'CHR' and 'POS' in tbl" = all(c("CHR", "POS") %in% colnames(tbl)))

  subset <- dplyr::slice_sample(tbl, n = {{ n_snps }})
  b38 <- map_to_dbsnp(tbl = subset, build = "38", by = "chr:pos", dbsnp_path = dbsnp_path)
  b37 <- map_to_dbsnp(tbl = subset, build = "37", by = "chr:pos", dbsnp_path = dbsnp_path)

  if(nrow(b37) > nrow(b38)) build <- "37" else build <- "38"
  cli::cli_inform("{nrow(b38)} snps matched GRCh38, {nrow(b37)} for GRCh37, inferring build to be {build}")

  build


}



# -------------------------------------------------------------------------


check_incompat_alleles <- function(tbl, dbsnp_df) {

  stopifnot(
    "Requires the following columns in tbl: rowid ,RSID, EffectAllele OtherAllele"  =
      all(c("rowid", "EffectAllele", "OtherAllele", "RSID") %in% colnames(tbl))
  )
  stopifnot(
    "Requires the following columns in dbsnp_df: CHR, POS, RSID, ref_allele,  alt_alleles"  =
      all(c("CHR","POS", "ref_allele", "alt_alleles", "RSID") %in% colnames(dbsnp_df))
  )


  # correct cols
  tbl <- dplyr::select(tbl, dplyr::all_of(c("rowid", "CHR", "POS", "RSID", "EffectAllele", "OtherAllele")))
  dbsnp_df <-
    dplyr::select(dbsnp_df, dplyr::all_of(c("CHR", "POS","RSID", "ref_allele", "alt_alleles"))) |>
    # multi-allelic SNPs into separate rows
    tidyr::separate_longer_delim(alt_alleles, delim =",")



  # check that EffectAllele ad OtherAllele is either REF ALT, or ALT REF
  no_strand_flip <- flag_incompat_alleles(tbl, dbsnp_df)
  failed <- dplyr::anti_join(tbl, no_strand_flip, by = "rowid")
  strand_flip <- flag_incompat_alleles(strand_flip(dplyr::anti_join(tbl, no_strand_flip, by = "rowid")), dbsnp_df)
  merged <- dplyr::bind_rows(no_strand_flip, strand_flip)


  # finished ----------------------------------------------------------------

  dplyr::mutate(tbl, incompat_alleles = dplyr::if_else(rowid %in% merged$rowid, FALSE, TRUE))


}


# -------------------------------------------------------------------------


flag_incompat_alleles <- function(tbl, dbsnp_df) {
  # if input contains
  dplyr::inner_join(tbl, dbsnp_df, by = c("CHR", "POS", "RSID"), relationship = "many-to-many") |>
    dplyr::mutate(
      a1_is_ref = dplyr::if_else(EffectAllele == ref_allele & OtherAllele == alt_alleles, TRUE, FALSE),
      a2_is_ref = dplyr::if_else(OtherAllele == ref_allele & EffectAllele == alt_alleles, TRUE, FALSE),
    ) |>
    dplyr::mutate(incompat_alleles = dplyr::if_else(a1_is_ref | a2_is_ref, FALSE, TRUE)) |>
    dplyr::select(-a1_is_ref, -a2_is_ref, -ref_allele, -alt_alleles) |>
    dplyr::filter(!incompat_alleles)
}


# -------------------------------------------------------------------------



flatten_dbsnp <- function(dbsnp_df) {

  stopifnot("flatten_dbsnp requires the following columns in dbsnp_df: CHR, POS, RSID, ref_allele,  alt_alleles"  =
              all(c("CHR","POS", "ref_allele", "alt_alleles", "RSID") %in% colnames(dbsnp_df))
  )
  # dbSNP contains IUPAC ambiguity codes for cases where the reference genome
  # can be multiple positions. (ex: ref_allele = "K", which means it can be a
  # a G or a T)
  # i want to map each ambiguity code (K, M, V etc, to a a vector of [A,C,G,T]),
  # separated by ","
  # so that i can easily expand the rows with tidyr::separate_longer_delim
  # for this, i need to convert the Biostrings::IUPAC_CODE_MAP
  # updated <-
  #   stringr::str_split(Biostrings::IUPAC_CODE_MAP, "") |>
  #   purrr::map(\(x) stringr::str_flatten(x, collapse=",")) |>
  #   purrr::map_chr(stringr::str_c) |>
  #   purrr::set_names(names(Biostrings::IUPAC_CODE_MAP))
  # use_data(updated, internal = TRUE)
  # stored in R/sysdata.rds


  # to deal with this, i expand the possible reference alleles to multiple rows,
  # allowing a SNP to match the specific ref allele.
  # handle the cases where the REF or ALT allele can be multiple allelles
  # by splitting the rows into separate rows. Each row will have one REF and one ALT.

  dbsnp_df |>
    dplyr::mutate(
      ref_allele = updated[ref_allele],
      alt_alleles = purrr::map_chr(alt_alleles, \(string) stringr::str_flatten(string, ","))
    ) |>
    tidyr::separate_longer_delim(alt_alleles, delim =",") |>
    tidyr::separate_longer_delim(ref_allele, delim =",")

}

# -------------------------------------------------------------------------





