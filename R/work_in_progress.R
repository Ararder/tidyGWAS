#
# if(!missing(dbsnp_path)) {
#
#   filter_funcs <-  purrr::map(paste0(filepaths$removed_rows,  c("main", "without_rsid"), "_incompat_alleles_or_not_in_dbsnp"), make_callback)
#   cli::cli_h2("Validating CHR/POS/RSID dbSNP for main rows")
#   main <- validate_with_dbsnp(data_list$main, build = build, dbsnp_path, filter_func = filter_funcs[[1]])
#
#   if(!is.null(data_list$without_rsid)) {
#
#     cli::cli_h2("Validating CHR/POS/RSID using dbSNP for rows without RSID")
#     data_list$without_rsid <- validate_with_dbsnp(data_list$without_rsid, build = build, dbsnp_path, filter_func = filter_funcs[[2]])
#
#
#     main <- dplyr::bind_rows(main,data_list$without_rsid)
#     before_unique_check <- dplyr::select(main, rowid)
#
#     # handle duplications
#     b38_missing <- dplyr::filter(main, is.na(CHR)) |>
#       dplyr::distinct(CHR_37,POS_37,EffectAllele,OtherAllele, .keep_all = TRUE)
#
#     main <- dplyr::filter(main, !is.na(CHR)) |>
#       dplyr::distinct(CHR,POS,EffectAllele,OtherAllele, .keep_all = TRUE) |>
#       dplyr::bind_rows(b38_missing)
#
#
#     removed <- dplyr::anti_join(before_unique_check, main, by = "rowid")
#     if(nrow(removed) > 0) {
#       cli::cli_alert_info("Found {nrow(removed)} rows which map to a CHR:POS:REF:ALT that another variant maps to. These are removed")
#       arrow::write_parquet(removed, paste0(filepaths$removed_rows, "map_to_dbsnp_duplications.parquet"))
#     }
#
#
#   }
# } else {
#   main <- data_list$main
# }
#
#
#
#
#
#
#
# validate_with_dbsnp <- function(data_list, build = c("NA", "37", "38"), dbsnp_path) {
#   rlang::check_required(dbsnp_path)
#   build <- rlang::arg_match(build)
#   data_list$main <- test_sumstat
#   tmp <- repair_dbnsp(data_list$main, dbsnp_path = dbsnp_path, build = build)
#
#   if(!is.null(data_list$without_rsid)) {
#     data_list$without_rsid) <-  repair_dbnsp(data_list$without_rsid, dbsnp_path = dbsnp_path, build = build)
#   }
#
#
#
#
#
# }
#
#
#
#
# repair_dbnsp <- function(tbl, dbsnp_path, build) {
#   # existence of chr:pos or rsid decides which columns to repair
#   has_rsid <- "RSID" %in% colnames(tbl)
#   has_chr_pos <- all(c("CHR", "POS") %in% colnames(tbl))
#
#   if(has_rsid & !has_chr_pos) {
#     cli::cli_inform("Using RSID to align with dbSNP")
#     tbl <- repair_chr_pos(tbl, dbsnp_path = dbsnp_path)
#
#   } else if(has_chr_pos & !has_rsid) {
#     cli::cli_inform("Using CHR and POS to align with dbSNP")
#     tbl <- repair_rsid(tbl, build = build, dbsnp_path =  dbsnp_path)
#
#
#   } else if(has_chr_pos & has_rsid) {
#     cli::cli_inform("Using CHR, POS and RSID to align with dbSNP")
#     tbl <- verify_chr_pos_rsid(tbl, build = build, dbsnp_path)
#
#   }
#
#   tbl
# }
