utils::globalVariables(c(
  "RSID", "POS", "CHR", "EffectAllele", "OtherAllele",
  "seqnames", "pos", "RefSNP_id", "ref_allele", "alt_alleles",
  "a1_is_ref", "a2_is_ref","incompat_alleles", "rowid", "RSID.x"
  ))

# https://www.ncbi.nlm.nih.gov/snp/docs/rs_multi_mapping/

# -------------------------------------------------------------------------
# This is the workhorse of performing joins between the ~1 billion rows of dbSNP
# and the inputted summary statistics. default argument for duplications = "smallest_rsid"
# which first arranges the data.frame on RSID, and keeps the first instance of each duplicate RSID


map_to_dbsnp <- function(tbl, build = c("37", "38"), by = c("rsid", "chr:pos"), duplications = c("smallest_rsid", "keep"), dbsnp_path) {


  # checks ------------------------------------------------------------------

  rlang::check_required(tbl)
  rlang::check_required(dbsnp_path)
  by = rlang::arg_match(by)
  build <- rlang::arg_match(build)
  duplications <- rlang::arg_match(duplications)

  if(!dplyr::is.tbl(tbl)) tbl <- dplyr::tibble(tbl)
  stopifnot("Cannot map 0 length tibble to dbsnp" = nrow(tbl) > 0)

  # open dbsnp dataset
  path_37 <- paste(dbsnp_path, "GRCh37", sep = "/")
  path_38 <- paste(dbsnp_path, "GRCh38", sep ="/")
  if(build == "37") path <- path_37 else path <- path_38
  dset <- arrow::open_dataset(path)


  # join with dbsnp ----------------------------------------------------------

  if(by == "rsid") {
    stopifnot("When mapping with RSID, RSID needs to present" = "RSID" %in% colnames(tbl))
    stopifnot(is.character(tbl$RSID))

    # RSID is stored as integer with "rs" suffix removed for better performance
    tbl$RSID <- as.integer(stringr::str_sub(tbl$RSID, start = 3))
    # batching by CHR gives SIGNIFICANTLY better performance
    # even though we check ALL rows against each chromosome
    chrom <- c(1:22, "X", "Y", "MT")

    res <-
      purrr::map(chrom, \(chrom) dplyr::filter(dset, CHR == {{ chrom }} & RSID %in% tbl$RSID) |> dplyr::collect()) |>
      purrr::list_rbind()

  } else {

    stopifnot("'CHR' and 'POS' need to be present in tbl" =  all(c("CHR", "POS") %in% colnames(tbl)))


    tbl$CHR <- as.character(tbl$CHR)
    tbl$POS <- as.integer(tbl$POS)

    res <-
      split(tbl, tbl$CHR) |>
      purrr::imap(\(df, chrom) dplyr::filter(dset, CHR == {{ chrom }} & POS %in% df$POS) |> dplyr::collect()) |>
      purrr::list_rbind()


  }



  dbsnp <- dplyr::rename(res, ref_allele = "REF", alt_alleles = "ALT") |>
    dplyr::mutate(RSID = stringr::str_c("rs", RSID)) |>
    # ensure correct types
    dplyr::mutate(CHR = as.character(CHR), POS = as.integer(POS), RSID = as.character(RSID)) |>
    # against dbSNP weirdness
    dplyr::distinct(CHR, POS, RSID, .keep_all = TRUE)


  # handle CHR:POS mapping to more than 1 RSID ------------------------------


  if(duplications == "smallest_rsid") {
    dbsnp <- dplyr::distinct(dplyr::arrange(dbsnp, RSID), CHR, POS, .keep_all = TRUE)
  }


  dbsnp

}





#' Use CHR and POS to get RSID from dbSNP v.155
#'
#' @description
#' [repair_rsid()] and [repair_chr_pos()] work similary, but assumes that either
#' rsid or chr_pos is missing.
#'
#' The functions first check which build the input data is in using [infer_build()].
#' Secondly, it maps each row to dbSNP using either rsid or chr_pos. Rows without a
#' match are flagged with `no_dbsnp_entry` = TRUE. Subsequently, rows where
#' EffectAllele and OtherAllele do not match the alleles in dbSNP are flagged with
#' `incompat_alleles` = TRUE.
#' Lastly, whichever build was missing is added to the data.
#'
#'
#' @param tbl a [dplyr::tibble()], formated with [tidyGWAS_columns()]
#' @inheritParams tidyGWAS
#' @return a [dplyr::tibble()] with columns `CHR`, `POS`, `POS_37` and `CHR_37` added
#'
#' @export
#'
#' @examples \dontrun{
#' sumstat_df <- repair_rsid(sumstat, bsgenome_list)
#' }
#'
repair_rsid <- function(tbl, build = c("NA", "37", "38"), dbsnp_path, add_missing_build=TRUE){

  # parse args -------------------------------------------------------------

  build = rlang::arg_match(build)
  if("RSID" %in% colnames(tbl)) {
    previous_rsid <- dplyr::select(tbl, c("rowid", "RSID"))
    tbl <- dplyr::select(tbl, -dplyr::all_of("RSID"))
  }
  if(build == "NA") build <- infer_build(tbl, dbsnp_path = dbsnp_path)



  # map to dbsnp ------------------------------------------------------------


  dbsnp_ref <- map_to_dbsnp(tbl, build = build, by = "chr:pos", dbsnp_path = dbsnp_path)

  # Add RSID, and flag any rows without a matching dbSNP entry
  tbl <-
    dplyr::left_join(tbl, dplyr::select(dbsnp_ref, CHR,POS, RSID, ref_allele), by = c("CHR","POS")) |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(RSID), TRUE, FALSE))


  # identify rows with incompatible alleles
  flags <- dplyr::select(check_incompat_alleles(dplyr::filter(tbl, !no_dbsnp_entry), dbsnp_ref), dplyr::all_of(c("rowid", "incompat_alleles")))

  # merge in flag
  tbl <- dplyr::left_join(tbl, flags, by = "rowid")

  # add missing build
  if(isTRUE(add_missing_build)) {
    tbl <- add_missing_build(tbl, ifelse(build == "38", "37", "38"), dbsnp_path = dbsnp_path)
  }

  tbl

}


#' Get CHR and POS using RSID
#' @inherit repair_rsid description params
#' @return a [dplyr::tibble()] with columns `CHR`, `POS`, `POS_37` and `CHR_37` added
#' @export
#'
#' @examples \dontrun{
#' sumstat_df <- repair_chr_pos(sumstat)
#' }
#'
repair_chr_pos <- function(tbl, dbsnp_path, add_missing_build = TRUE) {
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
#' @inheritParams repair_rsid
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





