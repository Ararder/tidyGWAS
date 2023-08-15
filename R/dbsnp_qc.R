utils::globalVariables(c(
  "RSID", "POS", "CHR", "EffectAllele", "OtherAllele",
  "seqnames", "pos", "RefSNP_id", "ref_allele", "alt_alleles",
  "a1_is_ref", "a2_is_ref","incompat_alleles", "rowid"
  ))

# https://www.ncbi.nlm.nih.gov/snp/docs/rs_multi_mapping/
# -------------------------------------------------------------------------

#' Compare CHR, POS and RSID with dbSNP reference data
#'
#' @param sumstat sumstat in tibble format, tidyGWAS column names
#' @param build Optional. Genome build, either 37 or 38. If not passed,
#' [infer_build()] will be called to get the genome build.
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' gwas <- tidyGWAS::test_file
#' bs <- get_bsgenome()
#' # make_callback can be passed a filepath to write out removed rows to.
#' callback <- make_callback("~/output_folder/verify_chr_pos_rsid_removed_rows.tsv")
#' verify_chr_pos_rsid(gwas, bs, build = 37)
#' }
verify_chr_pos_rsid <- function(sumstat, build = c("NA", "37", "38")) {

  # start -------------------------------------------------------------------
  build = rlang::arg_match(build)
  start_repair_message("verify_chr_pos_rsid")
  if(!"rowid" %in% colnames(sumstat)) sumstat$rowid <- 1:nrow(sumstat)
  if(build == "NA") build <- infer_build(sumstat)



  # get dbsnp data using CHR and POS
  dbsnp_ref <- map_to_dbsnp(sumstat, build = build, by = "chr:pos")
  # check for rows that doesn't match
  chr_pos_not_in_dbsnp <- dplyr::anti_join(sumstat, dbsnp_ref, by = c("CHR", "POS"))
  # update with RSID from dbsnp
  updated <- dplyr::inner_join(dplyr::select(sumstat, -RSID), dplyr::select(dbsnp_ref, -ref_allele, -alt_alleles), by = c("CHR", "POS"))
  # check how many rows had RSIDs different from what was found in dbSNP
  rsid_was_updated <- dplyr::inner_join(sumstat, dbsnp_ref, by = c("CHR", "POS")) |>
    dplyr::filter(RSID.x != RSID.y)

  # if any rows could not be found with CHR:POS mapping, try mapping with RSID
  if(nrow(chr_pos_not_in_dbsnp > 0)) {

    dbsnp_ref_2 <- map_to_dbsnp(chr_pos_not_in_dbsnp, build = build, by = "rsid")

    updated_2 <- dplyr::inner_join(
      dplyr::select(chr_pos_not_in_dbsnp, -CHR, -POS),
      dplyr::select(dbsnp_ref_2, -ref_allele, -alt_alleles),
      by = "RSID"
    )
  } else {
    updated_2 <- dplyr::tibble(RSID ="tmp", .rows = 0)
  }
  # merge the two dbsnp tables
  dbsnp_ref <- dplyr::bind_rows(dbsnp_ref, dbsnp_ref_2) |>
    dplyr::distinct()

  # merge the two rounds of mapping
  updated <- dplyr::bind_rows(updated, updated_2)

  rows_with_changed_values <- dplyr::bind_rows(rsid_was_updated, updated_2)
  removed <- dplyr::filter(sumstat, !rowid %in% updated$rowid)



  # -------------------------------------------------------------------------




  if(nrow(rows_with_changed_values) > 0) cli::cli_alert_info("A total of {nrow(rows_with_changed_values)} rows had incompatible CHR:POS and RSID data. These have been updated")


  # add flag that EffectAllele/OtherAllele is incompatible with REF/ALT in dbSNP
  cli::cli_li("Checking for incompatible alleles")
  flags <- check_incompat_alleles(updated, dbsnp_df = dbsnp_ref) |>
    dplyr::select(dplyr::all_of(c("rowid", "incompat_alleles")))



  sumstat <-
    dplyr::select(sumstat, rowid) |>
    dplyr::left_join(updated, by = "rowid") |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(CHR), TRUE, FALSE)) |>
    dplyr::left_join(flags, by = "rowid")


  #3) add CHR and POS from remaining build ------------------------------------
  add_missing_build(sumstat, ifelse(build == "37", "38", "37"))


}


#' Use CHR and POS to get RSID from dbSNP 155
#'
#'
#' @param sumstat a tibble with CHR,POS, EffectAllele and OtherAllele
#'
#' @param build Optional. Genome build, either NA, 37 or 38. If NA,
#' [infer_build()] will be called to get the genome build.
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' sumstat_df <- repair_rsid(sumstat, bsgenome_list)
#' }
#'
repair_rsid <- function(sumstat, build = c("NA", "37", "38")){
  build = rlang::arg_match(build)
  sumstat <- dplyr::select(sumstat, -any_of("RSID"))
  if(!"rowid" %in% colnames(sumstat)) sumstat$rowid <- 1:nrow(sumstat)
  if(build == "NA") build <- infer_build(sumstat)
  start_repair_message("repair_rsid")


  # Get RSID using CHR:POS
  dbsnp_ref <- map_to_dbsnp(sumstat, build = build, by = "chr:pos")

  # Add RSID, and flag any rows without a matching dbSNP entry
  sumstat <-
    dplyr::left_join(sumstat, dplyr::select(dbsnp_ref, CHR,POS, RSID), by = c("CHR","POS")) |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(RSID), TRUE, FALSE))

  # identify rows with incompatible alleles
  flags <- dplyr::select(check_incompat_alleles(sumstat, dbsnp_ref), dplyr::all_of(c("rowid", "incompat_alleles")))

  # merge in flag
  sumstat <- dplyr::left_join(sumstat, flags, by = "rowid")

  # add missing biuld
  add_missing_build(sumstat, ifelse(build == "38", "37", "38"))


}





#' Get CHR and POS using RSID
#'
#' @param sumstat a dplyr::tibble with with atleast RSID, EffectAllele and OtherAllele
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' sumstat_df <- repair_chr_pos(sumstat)
#' }
#'
repair_chr_pos <- function(sumstat) {
  sumstat <- dplyr::select(sumstat, -any_of(c("CHR", "POS")))
  if(!"rowid" %in% colnames(sumstat)) sumstat$rowid <- 1:nrow(sumstat)
  start_repair_message("repair_chr_pos")


  # start -------------------------------------------------------------------------

  # use RSID to get CHR and POS on GRCh38
  dbsnp_38 <- map_to_dbsnp(tbl = sumstat, build = "38", by = "rsid")

  # add CHR and POS, and flag that indicates if a match was found in dbsnp
  sumstat <- dplyr::left_join(sumstat, dplyr::select(dbsnp_38, CHR, POS, RSID), by = "RSID") |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(CHR), TRUE, FALSE))

  # add flags to see if REF/ALT is consistent with EA/OA
  flags <- check_incompat_alleles(sumstat, dbsnp_df = dbsnp_38) |>
    dplyr::select(rowid, incompat_alleles)


  # add GRCh37 CHR POS columns
  sumstat <- add_missing_build(sumstat, missing_build = "37")

  # return ----------------------------------------------------------------

  dplyr::left_join(sumstat, flags, by = "rowid")


}


# -------------------------------------------------------------------------


add_missing_build <- function(sumstat, missing_build = c("37", "38")) {
  missing_build = rlang::arg_match(missing_build)

  dbsnp <- map_to_dbsnp(sumstat, build = missing_build, by = "rsid") |>
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

make_callback <- function(outpath) {
  # update the filepath if it exists
  i <- 0
  while(file.exists(outpath)) outpath <- glue::glue(stringr::str_remove(outpath, "(_\\d{1})?\\.log.gz"), "_", (i <- i +1), ".log.gz")

  callback <- function(tbl) {
    # split into filter flags
    flags <- dplyr::select(tbl, rowid, dplyr::where(is.logical))
    if(ncol(flags) == 1) {
      cli::cli_inform("Found no flags to filter on")
      return(tbl)
    }

    remove <- dplyr::filter(flags, dplyr::if_any(dplyr::where(is.logical), \(x) x))
    count_by_flag <- purrr::map(dplyr::select(remove, -rowid), \(x) sum(x, na.rm = T))



    if(nrow(remove) > 0) {
      cli::cli_h2("Listing how many rows are removed per flag: ")
      cli::cli_dl(purrr::list_simplify(count_by_flag))
      cli::cli_inform("Removed a total of {nrow(remove)} rows: {.file {outpath}}")
    } else {
      cli::cli_h2("{.emph No rows were removed}")
    }



    data.table::fwrite(remove, outpath, sep = "\t")

    dplyr::select(tbl, -dplyr::where(is.logical)) |>
      dplyr::filter(!rowid %in% remove$rowid)
  }

  callback
}

# -------------------------------------------------------------------------


start_repair_message <-  function(func) {
  cli::cli_h1("tidyGWAS::{func}")
  cli::cli_ul()
  if(func == "repair_chr_pos")      cli::cli_li("Acquiring CHR and POS from dbSNP 155")
  if(func == "repair_rsid")         cli::cli_li("Acquiring RSID from dbSNP 155")
  if(func == "verify_chr_pos_rsid") cli::cli_li("Updating CHR:POS:RSID using reference data from dbSNP 155")

  cli::cli_li("{.code incompat_alleles} flags where EffectAllele/OtherAllele does not match REF/ALT")
  cli::cli_li("{.code no_dbsnp_entry} flags rows without dbSNP entry")

  cli::cli_li("This will likely take a few minutes...")

}









# -------------------------------------------------------------------------


infer_build <- function(sumstat, n_snps = 10000) {
  cli::cli_alert_info("Inferring build by checking {n_snps} snps matches against GRCh37 and GRCh38")
  stopifnot("Need 'CHR' and 'POS' in tbl" = all(c("CHR", "POS") %in% colnames(sumstat)))

  subset <- dplyr::slice_sample(sumstat, n = {{ n_snps }})
  b38 <- map_to_dbsnp(tbl = subset, build = "38", by = "chr:pos")
  b37 <- map_to_dbsnp(tbl = subset, build = "37", by = "chr:pos")

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
  tbl <- dplyr::select(tbl, dplyr::all_of(c("rowid", "RSID", "EffectAllele", "OtherAllele")))
  dbsnp_df <-
    dplyr::select(dbsnp_df, dplyr::all_of(c("RSID", "ref_allele", "alt_alleles"))) |>
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
  dplyr::inner_join(tbl, dbsnp_df, by = "RSID") |>
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





