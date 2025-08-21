utils::globalVariables(c("CHR_38","ALT"))
map_indels_dbsnp <- function(tbl, by = "chr:pos", dbsnp_path) {

  if(by == "rsid") {
    id  <- "RSID"
    tbl <- dplyr::select(tbl, -dplyr::any_of(c("CHR", "POS")))

    # aligning to both builds
    p1 <- fs::path(dbsnp_path, "dbSNP157_indels", "37")
    p2 <- fs::path(dbsnp_path, "dbSNP157_indels", "38")
    cli::cli_inform("looking for dbsnp files in {p1} and {p2}")
    file.exists(p1) || cli::cli_abort("Could not find the expected indel reference data. Have you updated the reference data?")
    file.exists(p2) || cli::cli_abort("Could not find the expected indel reference data. Have you updated the reference data?")
    dset_37 <- arrow::open_dataset(p1)
    dset_38 <- arrow::open_dataset(p2)

    ref_37 <- match_by(tbl, "RSID", dset_37)
    ref_38 <- match_by(tbl, "RSID", dset_38)



    ref <- dplyr::inner_join(ref_37, ref_38, by = c("RSID", "EffectAllele", "OtherAllele"), suffix = c("_37", "_38")) |>
      dplyr::filter(!is.na(CHR_37) & !is.na(CHR_38)) |>
      dplyr::filter(CHR_37 == CHR_38) |>
      dplyr::rename(CHR = CHR_38) |>
      dplyr::select(-dplyr::any_of(c("CHR_37")))

    merged <- dplyr::inner_join(tbl, ref, by = c("RSID", "EffectAllele", "OtherAllele"))



  } else {

    tbl <- dplyr::select(tbl, -dplyr::any_of(c("RSID")))
    build <- get_build(tbl, dbsnp_path)
    dset <- arrow::open_dataset(fs::path(dbsnp_path, "dbSNP157_indels", build))

    ref_data <- match_by(tbl, c("CHR", "POS"), dset)


    other_build <- ifelse(build == "38", "37", "38")
    other_path <- arrow::open_dataset(fs::path(dbsnp_path, "dbSNP157_indels",  other_build))
    other_b <- match_by(dplyr::select(ref_data, RSID, OtherAllele, EffectAllele), "RSID", other_path)

    suffix  <- c(paste0("_", build), paste0("_", other_build))

    # in both builds
    ref  <- dplyr::inner_join(ref_data, other_b, by = c("RSID", "EffectAllele", "OtherAllele"), suffix = suffix)

    merged <-
      dplyr::inner_join(tbl, ref, by = c("CHR" = paste0("CHR_", build), "POS" = paste0("POS_", build), "EffectAllele", "OtherAllele")) |>
      dplyr::select(-dplyr::any_of(c(paste0("CHR_", other_build)))) |>
      dplyr::rename_with(\(x) paste0("POS_", build), "POS")


  }

  cli::cli_alert_info(
    "Removed a total of {nrow(tbl) - nrow(merged)} variants.
		Possible reasons for removal:
		(1) Could not find a match in dbSNP
		(2) Different chromosomes across builds
		(3) Not present in both builds"
  )

  merged

}




match_by <- function(tbl, vec, dset) {

  vec <- c(vec, "REF", "ALT")

  ref_ALT_EFF <- dset |>
    dplyr::semi_join(
      dplyr::rename(tbl, REF = OtherAllele, ALT = EffectAllele),
      by = vec
    ) |>
    dplyr::collect()


  ref_REF_EFF <- dset |>
    dplyr::semi_join(
      dplyr::rename(tbl, REF=EffectAllele, ALT = OtherAllele),
      by = vec
    ) |>
    dplyr::collect()


  cli::cli_inform("Found {nrow(ref_ALT_EFF)} matches using ALT as EffectAllele, and {nrow(ref_REF_EFF)} matches using REF as EffectAllele")
  if(nrow(ref_ALT_EFF) > nrow(ref_REF_EFF)) {
    dplyr::mutate(ref_ALT_EFF, tmp = REF) |>
      dplyr::rename(EffectAllele = ALT, OtherAllele = REF) |>
      dplyr::rename(REF = tmp)

  } else {
    dplyr::mutate(ref_REF_EFF, tmp = REF) |>
      dplyr::rename(EffectAllele = REF, OtherAllele = ALT) |>
      dplyr::rename(REF = tmp)

  }

}

get_build <- function(tbl, dbsnp_path, n_snps = 10000) {
  subset <- dplyr::slice_sample(dplyr::filter(tbl, CHR == "21"), n = {{ n_snps }})
  if(nrow(subset) != n_snps) subset <- dplyr::slice_sample(tbl, n = {{ n_snps }})


  p37 <- fs::path(dbsnp_path, "dbSNP157_indels","37")
  p38 <- fs::path(dbsnp_path, "dbSNP157_indels","38")
  file.exists(p37) || cli::cli_abort("Could not find the expected indel reference data. Have you updated the reference data?")
  file.exists(p38) || cli::cli_abort("Could not find the expected indel reference data. Have you updated the reference data?")

  b37 <- match_by(subset, c("CHR", "POS"), arrow::open_dataset(p37))
  b38 <- match_by(subset, c("CHR", "POS"), arrow::open_dataset(p38))

  if(nrow(b38) > nrow(b37)) {
    cli::cli_alert_info("Build detected as GRCh38")
    return("38")
  } else {
    cli::cli_alert_info("Build detected as GRCh37")
    return("37")
  }


}
