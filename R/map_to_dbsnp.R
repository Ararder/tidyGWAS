# map_to_dbsnp contains the code to query dbSNP, either using RSID, or chromosome ("CHR")
# and position ("POS"). tidyGWAS supports two ways of interacting with dbSNP:
# 1) through the "arrow" implementation, which has hive-style partitioned
# dbSNP v155 in parquet gzipped files
# 2) through the BSgenome package from bioconductor, which represents dbSNP v155
# in their own internal format

map_to_dbsnp <- function(tbl, build = c("37", "38"), by = c("rsid", "chr:pos"), implementation = c("arrow", "bsgenome"), duplications = c("smallest_rsid", "keep")) {
  by = rlang::arg_match(by)
  implementation = rlang::arg_match(implementation)
  build <- rlang::arg_match(build)
  duplications <- rlang::arg_match(duplications)


  map_to_dbsnp_checks(tbl = tbl, build = build, by = by, implementation = implementation)


  # arrow or bsgenome ----------------------------------------------------

  if(implementation == "arrow") {
    dbsnp <- map_to_dbsnp_arrow(tbl = tbl, build = build, by= by)
  } else {
    dbsnp <- map_to_dbsnp_bsgenome(tbl = tbl, build = build, by= by)
  }


  # ensure correct types
  tmp <- dplyr::mutate(dbsnp, CHR = as.character(CHR), POS = as.integer(POS), RSID = as.character(RSID)) |>
    dplyr::distinct(CHR, POS, RSID, .keep_all = TRUE)

  if(duplications == "smallest_rsid") {
    tmp <- dplyr::distinct(dplyr::arrange(tmp, RSID), CHR, POS, .keep_all = TRUE)
  }


  tmp

}

map_to_dbsnp_checks <- function(tbl, build, by, implementation = "arrow") {


  stopifnot("tbl should be a tibble" = "tbl" %in% class(tbl))
  stopifnot("tbl is empty" = nrow(tbl) > 0)
  if(by == "rsid") {
    stopifnot("RSID" %in% colnames(tbl))
    stopifnot(is.character(tbl$RSID))
  } else {
    stopifnot("'CHR' and 'POS' need to be present in tbl" =  all(c("CHR", "POS") %in% colnames(tbl)))
  }





}


map_to_dbsnp_arrow <- function(tbl, build, by) {

  if(build == "37") path <- Sys.getenv("grch37") else path <- Sys.getenv("grch38")
  dset <- arrow::open_dataset(path)


  if(by == "chr:pos") {

    tbl$CHR <- as.character(tbl$CHR)
    tbl$POS <- as.integer(tbl$POS)

    res <-
      split(tbl, tbl$CHR) |>
      purrr::imap(\(df, chrom) dplyr::filter(dset, CHR == {{ chrom }} & POS %in% df$POS) |> dplyr::collect()) |>
      purrr::list_rbind()


  } else if(by == "rsid") {

    tbl$RSID <- as.integer(stringr::str_sub(tbl$RSID, start = 3))
    chrom <- c(1:22, "X", "Y", "MT")
    # batching by CHR gives SIGNIFICANTLY better performance
    # even though we run all rows per chrom in each %in%

    res <-
      c(1:22, "X", "Y", "MT") |>
      purrr::map(\(chrom) dplyr::filter(dset, CHR == {{ chrom }} & RSID %in% tbl$RSID) |> dplyr::collect()) |>
      purrr::list_rbind()


  }

  dplyr::rename(res, ref_allele = "REF", alt_alleles = "ALT") |>
    dplyr::mutate(RSID = stringr::str_c("rs", RSID))

}

map_to_dbsnp_bsgenome <- function(tbl, build, by) {



  bsgenome_objects <- get_bsgenome()
  if(build == "37") {
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
    dplyr::select(-dplyr::any_of(c("strand", "alleles_as_ambig", "genome_compat")))
    # this will ensure correct types AND that these columns exist
    # while possible superfluous, ensures type safety from output of BSgenome::snpsByXX functions




    dbsnp
}


get_bsgenome <- function() {
  list(
    "snps_37" = SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37,
    "genome_37" = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    "snps_38" = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
    "genome_38" = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
  )

}
