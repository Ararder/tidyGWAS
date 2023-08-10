bsgenome_checks <- function(tbl, build, by, implementation = "bsgenome") {

  stopifnot("Can only map using 'rsid' or 'chr:pos'" = by %in% c("rsid", "chr:pos"))
  stopifnot("implementation can only be arrow or bsgenome" = implementation %in% c("arrow", "bsgenome"))
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




}

map_to_dbsnp <- function(tbl, build = 37, by = "rsid", implementation = "bsgenome") {

  bsgenome_checks(tbl = tbl, build = build, by = by)


  # arrow or bsgenome impl ----------------------------------------------------

  if(implementation == "arrow") {
    dbsnp <- map_to_dbsnp_arrow(tbl = tbl, build = build, by= by)
  } else {
    dbsnp <- map_to_dbsnp_bsgenome(tbl = tbl, build = build, by= by)
  }


  # ensure correct types
  dplyr::mutate(dbsnp, CHR = as.character(CHR), POS = as.integer(POS), RSID = as.character(RSID)) |>
    dplyr::distinct(CHR, POS, RSID, .keep_all = TRUE)

}

map_to_dbsnp_arrow <- function(tbl, build, by, path_b37, path_b38) {

  b37 <- arrow::open_dataset("/nas/depts/007/sullilab/shared/arvhar/snp_level_annotatations/dnSNP155/GRCh37")
  b38 <- arrow::open_dataset("/nas/depts/007/sullilab/shared/arvhar/snp_level_annotatations/dnSNP155/GRCh38")

  if(by = "rsid") key <- "RSID" else key <- c("CHR", "POS")

  tbl$RSID <- as.integer(stringr::str_sub(tbl$RSID, start = 3))
  tbl$CHR <- as.character(tbl$CHR)
  tbl$POS <- as.integer(tbl$POS)
  out <- dplyr::semi_join(b37, tbl, by = key)

}

map_to_dbsnp_bsgenome <- function(tbl, build, by) {



  bsgenome_objects <- get_bsgenome()
  if(build == 37) {
    snps <- bsgenome_objects$snps_37
    genome <- bsgenome_objects$genome_37
  } else {
    snps <- bsgenome_objects$snps_38
    genome <- bsgenome_objects$genome_38
  }



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




    dbsnp
}

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
