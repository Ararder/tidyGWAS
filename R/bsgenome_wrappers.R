# Todo - what happens if input a RSID without a hit?
utils::globalVariables(c(
  "RSID", "POS", "CHR", "EffectAllele", "OtherAllele",
  "seqnames", "pos", "RefSNP_id", "ref_allele", "alt_alleles",
  "a1_is_ref", "a2_is_ref","uncompatible_alleles", "rowid"
  ))

# -------------------------------------------------------------------------
#' Check that CHR POS and RSID all make sense in a GWAS summary statistic
#'
#' @param sumstat sumstat in tibble format, tidyGWAS column names
#' @param bsgenome_objects use get_bsgenome()
#' @param build optional, can be used to skip the infer_build step
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
verify_chr_pos_rsid <- function(sumstat, bsgenome_objects, build) {

  # start -------------------------------------------------------------------

  start_repair_message("verify_chr_pos_rsid")
  if(!"rowid" %in% colnames(sumstat)) sumstat$rowid <- 1:nrow(sumstat)
  dbsnp <- vector("list")


  # 1) figure out genome build -------------------------------------------------

  if(missing(build)) build <- infer_build(sumstat, bsgenome_objects = bsgenome_objects)
  if(missing(bsgenome_objects)) bsgenome_objects <- get_bsgenome()

  # check if RSIDs agree
  dbsnp[[as.character(build)]] <- map_to_dbsnp(dplyr::select(sumstat, -RSID), build = build, by = "chr:pos", bsgenome_objects)

  # there are multiple cases where the same CHR:POS maps to multiple RSIDs
  # in that case, we select the first RSID (smallest rs number, since rsid is arranged by )
  reduced <- dplyr::distinct(dplyr::arrange(dbsnp[[as.character(build)]], RSID), CHR, POS, .keep_all = TRUE)

  # find the cases where chr_pos could not identify a RSID
  no_match <- dplyr::anti_join(sumstat, reduced, by = c("CHR", "POS"))
  first_round_match <- dplyr::inner_join(dplyr::select(sumstat, -RSID), reduced, by = c("CHR", "POS"))

  # check if the provided RSID exists in dbSNP
  if(nrow(no_match > 0)) {
    dbsnp[["rsid_mapping"]] <- map_to_dbsnp(dplyr::select(no_match, -CHR, -POS), build = build, by = "rsid", bsgenome_objects)
    reduced_round2 <- dplyr::distinct(dplyr::arrange(dbsnp[["rsid_mapping"]], RSID), CHR, POS, .keep_all = TRUE)
    second_round_match <- dplyr::inner_join(dplyr::select(no_match, -CHR, -POS), reduced_round2, by = "RSID")

  } else {
    second_round_match <- dplyr::tibble(RSID ="tmp", .rows = 0)
  }

  # collect

  merged <- dplyr::bind_rows(first_round_match, second_round_match)
  merged <- dplyr::distinct(merged, CHR, POS, EffectAllele, OtherAllele, .keep_all = TRUE)


  # find removed rows and updated rows
  removed <- dplyr::filter(sumstat, !rowid %in% merged$rowid)
  updated_rows <- dplyr::filter(sumstat, rowid %in% merged$rowid) |>
    dplyr::anti_join(merged, by = c("rowid", "CHR", "POS","RSID"))

  if(nrow(updated_rows) > 0) cli::cli_alert_info("A total of {nrow(updated_rows)} rows had incompatible CHR:POS and RSID data. These have been updated")


  # add flag that EffectAllele/OtherAllele is incompatible with REF/ALT in dbSNP
  flags <- qc_with_dbsnp(
    dplyr::select(merged, -ref_allele, -alt_alleles),
    dbsnp_df = dplyr::select(merged, -EffectAllele, -OtherAllele,rowid)
  ) |>
    dplyr::select(dplyr::all_of(c("rowid", "uncompatible_alleles")))


  # clean up
  sumstat <- dplyr::select(sumstat, rowid) |>
    dplyr::left_join(merged, by = "rowid") |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(CHR), TRUE, FALSE)) |>
    dplyr::left_join(flags, by = "rowid")



  #3) add CHR and POS from remaining build ------------------------------------


  add_build <- ifelse(build == 38, 37, 38)
  cli::cli_alert_info("CHR and POS were on GRCh{build}. Acquiring positions on GRCh{add_build}, by mapping to dbSNP with RSID")
  dbsnp[[as.character(add_build)]] <- map_to_dbsnp(dplyr::filter(sumstat, !is.na(RSID)), build = add_build, by = "rsid", bsgenome_objects = bsgenome_objects) |>
    dplyr::select(-dplyr::all_of(c("ref_allele", "alt_alleles"))) |>
    # remove multi-allelic SNPs
    dplyr::distinct(CHR, POS, RSID)


  if(add_build == 38) {
    final <- dplyr::left_join(dplyr::rename(sumstat, CHR_37 = CHR, POS_37 = POS), dbsnp[[as.character(add_build)]], by = "RSID")
  } else {
    final <- dplyr::left_join(sumstat, dplyr::rename(dbsnp[[as.character(add_build)]], CHR_37 = CHR, POS_37 = POS), by = "RSID")
  }



  # return ------------------------------------------------------------------


  dplyr::select(final, -dplyr::any_of(c("ref_allele", "alt_alleles")))

}


#' Get RSID from either GRCh37 or GRCh38 reference genome, using CHR and POS
#'
#' @param sumstat a dplyr::
#' tibble with atleast CHR,POS, EffectAllele and OtherAllele
#' @param bsgenome_objects a list containing BSgenome genoms and snp_locs. see get_bsgenome_objects
#' @param build genome build, either 37 or 38
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' sumstat_df <- repair_rsid(sumstat, bsgenome_list)
#' }
#'
repair_rsid <- function(sumstat, bsgenome_objects, build){
  if(!"rowid" %in% colnames(sumstat)) sumstat$rowid <- 1:nrow(sumstat)
  start_repair_message("repair_rsid")
  dbsnp <- vector("list")
  if(missing(bsgenome_objects)) bsgenome_objects <- get_bsgenome()


  # 1) figure out genome build -------------------------------------------------
  if(missing(build)) build <- infer_build(sumstat, bsgenome_objects=bsgenome_objects)

  # 2) get RSID and flags  --------------------------------------------------------
  dbsnp[[as.character(build)]] <- map_to_dbsnp(sumstat, build = build, by = "chr:pos", bsgenome_objects)

  # CHR:POS can map to multiple SNPs, select the
  tmp_rsid <- dplyr::distinct(dplyr::arrange(dbsnp[[as.character(build)]], RSID), CHR, POS, .keep_all = TRUE)

  # merge in RSID, and add flag any rows that could not find a RSID
  sumstat <- dplyr::left_join(sumstat, tmp_rsid, by = c("CHR","POS")) |>
    dplyr::mutate(no_dbsnp_entry = dplyr::if_else(is.na(RSID), TRUE, FALSE))

  # add flag that EffectAllele/OtherAllele is incompatible with REF/ALT in dbSNP
  flags <- qc_with_dbsnp(tbl = dplyr::filter(sumstat, !no_dbsnp_entry), dbsnp_df = dbsnp[[as.character(build)]]) |>
    dplyr::select(dplyr::all_of(c("rowid", "uncompatible_alleles")))

  # merge in flag
  sumstat <- dplyr::left_join(sumstat, flags, by = "rowid")



  #3) add CHR and POS from remaining build ------------------------------------


  add_build <- ifelse(build == 38, 37, 38)
  cli::cli_alert_info("CHR and POS were on GRCh{build}. Acquiring positions on GRCh{add_build}, by mapping to dbSNP with RSID")
  dbsnp[[as.character(add_build)]] <- map_to_dbsnp(sumstat, build = add_build, by = "rsid", bsgenome_objects = bsgenome_objects) |>
    dplyr::select(-dplyr::all_of(c("ref_allele", "alt_alleles"))) |>
    # remove multi-allelic SNPs
    dplyr::distinct(CHR, POS, RSID)


  if(add_build == 38) {
    final <- dplyr::left_join(dplyr::rename(sumstat, CHR_37 = CHR, POS_37 = POS), dbsnp[[as.character(add_build)]], by = "RSID")
  } else {
    final <- dplyr::left_join(sumstat, dplyr::rename(dbsnp[[as.character(add_build)]], CHR_37 = CHR, POS_37 = POS), by = "RSID")
  }


  # return --------------------------------------------------------------------


  dplyr::select(final, -dplyr::any_of(c("ref_allele", "alt_alleles")))

}


# -------------------------------------------------------------------------



#' Get CHR and POS from a reference genome with rsID
#'
#' @param sumstat a dplyr::tibble with with atleast RSID, EffectAllele and OtherAllele
#' @param bsgenome_objects a list containing BSgenome genoms and snp_locs. see get_bsgenome_objects
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' sumstat_df <- repair_chr_pos(sumstat, bsgenome_list)
#' }
#'
repair_chr_pos <- function(sumstat, bsgenome_objects){

  if(!"rowid" %in% colnames(sumstat)) sumstat$rowid <- 1:nrow(sumstat)
  start_repair_message("repair_chr_pos")
  if(missing(bsgenome_objects)) bsgenome_objects <- get_bsgenome()


  # use RSID to get CHR and POS on GRCh38
  dbsnp_38 <- map_to_dbsnp(sumstat, build = 38, by = "rsid", bsgenome_objects =  bsgenome_objects)

  # find duplicate
  dup_vec <- !(duplicated(dbsnp_38[,1:2]) | duplicated(dbsnp_38[,1:2], fromLast = TRUE))
  if(sum(!dup_vec) > 0) cli::cli_alert_warning("Found {sum(!dup_vec)} rows where the same CHR:POS maps to different RSIDs")


  # add flag to indicate whether RSID was found in dbSNP
  sumstat <- dplyr::mutate(sumstat, no_dbsnp_entry = dplyr::if_else(!RSID %in% dbsnp_38$RSID, TRUE, FALSE))

  # add flags to see if REF/ALT is consistent with EA/OA
  flags <- qc_with_dbsnp(tbl = dplyr::filter(sumstat, no_dbsnp_entry), dbsnp_df = dbsnp_38) |>
    dplyr::select(dplyr::any_of(c("rowid", "RSID", "uncompatible_alleles")))


  # get CHR and POS from GRCh37 ---------------------------------------------
  dbsnp_37 <- map_to_dbsnp(sumstat, build = 37, by = "rsid", bsgenome_objects = bsgenome_objects) |>
    dplyr::select(CHR_37 = CHR, POS_37 = POS, RSID) |>
    dplyr::distinct()

  # some weird duplications. Make sure none happen when adding CHR, POS,RSID
  dbsnp_38 <- dplyr::select(dbsnp_38, CHR, POS, RSID) |>
    dplyr::distinct()


  # merge with GRCh37 CHR and POS
  final <-
    # POS and CHR from GRCh38
    dplyr::left_join(sumstat, dplyr::select(dbsnp_38, CHR, POS, RSID), by = "RSID") |>
    # POS and CHR from GRCh37
    dplyr::left_join(dbsnp_37, by = "RSID") |>
    dplyr::left_join(flags, by = c("rowid", "RSID"))



  # return ----------------------------------------------------------------

  final

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
  if(func == "repair_chr_pos") {
    cli::cli_li("{.pkg BSgenome} package is used to acquire CHR and POS, with RSID as input")

  } else if(func == "repair_rsid") {
    cli::cli_li("{.pkg BSgenome} package is used to acquire RSID, with CHR and POS as input")

  } else if(func == "verify_chr_pos_rsid") {
    cli::cli_li("Uses CHR and POS to get RSID. For rows where CHR:POS could not find RSID,RSID is used to get CHR:POS")
  }
  cli::cli_li("{.code incompatible_alleles} flags where EffectAllele/OtherAllele does not match REF/ALT")
  cli::cli_li("{.code no_dbsnp_entry} flags rows without dbSNP entry")

  cli::cli_li("This will likely take a few minutes...")

}


# -------------------------------------------------------------------------


map_to_dbsnp <- function(tbl, build = 37, by = "rsid", bsgenome_objects) {

  # validate input types ----------------------------------------------------


  stopifnot("Can only map using 'rsid' or 'chr:pos'" = by %in% c("rsid", "chr:pos"))
  stopifnot("tbl should be a tibble" = "tbl" %in% class(tbl))
  stopifnot("tbl is empty" = nrow(tbl) > 0)
  if(by == "rsid") {
    stopifnot("RSID" %in% colnames(tbl))
    stopifnot(is.character(tbl$RSID))

  } else {
    stopifnot("'CHR' and 'POS' need to be present in tbl" =  all(c("CHR", "POS") %in% colnames(tbl)))
    stopifnot(is.character(tbl$CHR) & is.integer(tbl$POS))
  }
  stopifnot("Only supports GRCh37 or GRCh38" = build %in% c(37, 38))



  # attempt to use preloaded bsgenome ---------------------------------------


  if(missing(bsgenome_objects)) bsgenome_objects <- get_bsgenome()

  if(build == 37) {
    snps <- bsgenome_objects$snps_37
    genome <- bsgenome_objects$genome_37
  } else {
    snps <- bsgenome_objects$snps_38
    genome <- bsgenome_objects$genome_38
  }




  # join with dbsnp ---------------------------------------------------------

  if(by == "rsid") {

    res <- BSgenome::snpsById(
      x = snps,
      genome = genome,
      # make sure only valid RSIDs are passed to snpsById
      ids = dplyr::filter(tbl, stringr::str_detect(RSID, "rs\\d{1,9}"))[["RSID"]],
      ifnotfound="drop"
    ) |>
      data.table::as.data.table()

  } else {

    res <-
      # splitting by chromosome reduces memory usage of BSgenome::snpsByOverlaps
      split(tbl, tbl$CHR) |>
      purrr::map(\(chr_df) dplyr::arrange(chr_df, POS)) |>
      # convert to GPos object
      purrr::map(\(chr_df) GenomicRanges::GPos(chr_df[["CHR"]], chr_df[["POS"]])) |>
      purrr::map(\(chr_df_as_gpos) BSgenome::snpsByOverlaps(ranges = chr_df_as_gpos, x = snps, genome = genome)) |>
      purrr::map(data.table::as.data.table) |>
      # remove empty entries
      purrr::keep(\(x) nrow(x) > 0) |>
      purrr::list_rbind()
    if(rlang::is_empty(res)) res <- data.table::data.table("seqnames" = character(), "pos" = integer(), "RefSNP_id" = character(), "ref_allele" = character(),"alt_alleles" = list())

  }


  # clean up and return results ---------------------------------------------


  dbsnp <- res |>
    tibble::as_tibble() |>
    dplyr::rename(CHR = seqnames, POS = pos, RSID = RefSNP_id) |>
    dplyr::select(-dplyr::any_of(c("strand", "alleles_as_ambig", "genome_compat"))) |>
    # this will ensure correct types AND that these columns exist
    # while possible superfluous, ensures type safety from output of BSgenome::snpsByXX functions
    dplyr::mutate(
      CHR = as.character(CHR),
      POS = as.integer(POS),
      RSID = as.character(RSID)
    ) |>
    dplyr::distinct(CHR, POS, RSID, .keep_all = TRUE)



  dbsnp

}


# -------------------------------------------------------------------------


infer_build <- function(sumstat, n_snps = 10000,  ...) {
  cli::cli_alert_info("Inferring build by checking {n_snps} snps matches against GRCh37 and GRCh38")
  stopifnot("Need 'CHR' and 'POS' in tbl" = all(c("CHR", "POS") %in% colnames(sumstat)))

  subset <- dplyr::slice_sample(sumstat, n = {{ n_snps }})
  b38 <- map_to_dbsnp(tbl = subset, build = 38, by = "chr:pos", ...)
  b37 <- map_to_dbsnp(tbl = subset, build = 37, by = "chr:pos", ...)

  if(nrow(b37) > nrow(b38)) build <- 37 else build <- 38
  cli::cli_inform("{nrow(b38)} snps matched GRCh38, {nrow(b37)} for GRCh37, inferring build to be {build}")

  build


}



# -------------------------------------------------------------------------


qc_with_dbsnp <- function(tbl, dbsnp_df) {

  stopifnot("Requires the following columns in tbl: rowid ,RSID, EffectAllele OtherAllele"  =
              all(c("rowid", "EffectAllele", "OtherAllele", "RSID") %in% colnames(tbl))
  )
  stopifnot("Requires the following columns in dbsnp_df: CHR, POS, RSID, ref_allele,  alt_alleles"  =
              all(c("CHR","POS", "ref_allele", "alt_alleles", "RSID") %in% colnames(dbsnp_df))
  )


  # make sure input format is in correct format
  tbl <-      dplyr::select(tbl,      dplyr::all_of(c("rowid", "RSID", "EffectAllele", "OtherAllele")))
  dbsnp_df <- dplyr::select(dbsnp_df, dplyr::all_of(c("RSID", "CHR", "POS", "ref_allele", "alt_alleles")))
  # separate multi-allelic SNPs into separate rows
  dbsnp_df <- flatten_dbsnp(dbsnp_df)



  # first pass QC -----------------------------------------------------------
  # Apply no strand flipping, ignore ambigious SNPs. Try the naive way,
  # match by RSID and check that alleles make sense

  first_pass <-
    dplyr::left_join(tbl, dbsnp_df, by = "RSID") |>
    # check if ref allele corresponds to EffectAllele
    # also checks that EA/OA pairing is consistent with REF/ALT in dbsnp
    dplyr::mutate(
      a1_is_ref = dplyr::if_else(EffectAllele == ref_allele & OtherAllele == alt_alleles, TRUE, FALSE),
      a2_is_ref = dplyr::if_else(OtherAllele == ref_allele & EffectAllele == alt_alleles, TRUE, FALSE),
    ) |>
    # create a flag for inconsistency between EA/OA and ref/Alt
    # since multi-allelic SNPs are split, many of them are flagged as inconsistent
    dplyr::mutate(uncompatible_alleles = dplyr::if_else(a1_is_ref | a2_is_ref, TRUE, FALSE))

  # create a subset listing all rows with EA & OA that matches REF/ALT
  a1_a2_in_ref_or_alt <-
    dplyr::filter(first_pass, uncompatible_alleles) |>
    dplyr::select(rowid, RSID, a1_is_ref, uncompatible_alleles) |>
    dplyr::mutate(strand_flip = FALSE) |>
    # there are some cases where a multi-allelic SNP matches, and can be both
    # a1_is_ref, and a2_is_ref, see rs2765306 for example
    # remove this by calling distinct
    dplyr::distinct(rowid, RSID, .keep_all = TRUE)

  removed <-
    dplyr::anti_join(tbl, a1_a2_in_ref_or_alt, by = "RSID") |>
    dplyr::mutate(reason = "A1_A2 not in REF/ALT")

  # second pass qc ----------------------------------------------------------
  # It is possible that some SNPs are removed because they are on the wrong strand
  # they need to be strand flipped to align with the reference genome

  a1_a2_in_ref_or_alt_strand_flipped <-
    # strand flip and remove reason column
    dplyr::select(strand_flip(removed), rowid, RSID, EffectAllele, OtherAllele) |>
    dplyr::left_join(dbsnp_df, by = "RSID") |>
    dplyr::mutate(
      a1_is_ref = dplyr::if_else(EffectAllele == ref_allele & OtherAllele == alt_alleles, TRUE, FALSE),
      a2_is_ref = dplyr::if_else(OtherAllele == ref_allele & EffectAllele == alt_alleles, TRUE, FALSE),
    ) |>
    # Here any rows that don't match ref/ALT in dbSNP is removed
    dplyr::mutate(uncompatible_alleles = dplyr::if_else(a1_is_ref | a2_is_ref, TRUE, FALSE)) |>
    dplyr::select(rowid,RSID, a1_is_ref, uncompatible_alleles) |>
    dplyr::mutate(strand_flip = TRUE) |>
    dplyr::filter(uncompatible_alleles)

  # add flags to input tbl --------------------------------------------------
  out <- dplyr::left_join(
    tbl,
    dplyr::bind_rows(a1_a2_in_ref_or_alt, a1_a2_in_ref_or_alt_strand_flipped),
    by = c("RSID", "rowid")
  )

  stopifnot("more rows after qc than before" = nrow(out) == nrow(tbl))
  dplyr::mutate(out, uncompatible_alleles = !uncompatible_alleles)

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
  updated <-
    stringr::str_split(Biostrings::IUPAC_CODE_MAP, "") |>
    purrr::map(\(x) stringr::str_flatten(x, collapse=",")) |>
    purrr::map_chr(stringr::str_c) |>
    purrr::set_names(names(Biostrings::IUPAC_CODE_MAP))


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



#' Read in RsMergeArch file
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' rs_merge_arch <- get_ref_data()
#' }
get_ref_data <- function() {

  data.table::fread(Sys.getenv("rs_merge_arch")) |> dplyr::tibble()
}


# -------------------------------------------------------------------------


#' load SNPlocs and REF genome for GRCh 37 and 38
#'
#' @return a list of SNPlocs and BSgenome objects
#' @export
#'
#' @examples \dontrun{
#' bsgenome <- get_bsgenome()
#' }
get_bsgenome <- function() {
  list(
    SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37,
    BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
    BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    # BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  ) |>
    purrr::set_names(c("snps_37", "genome_37", "snps_38", "genome_38"))

}
