utils::globalVariables(c(
  "chr_pos", "n", "n_chr_pos",
  ".data", "invalid_rsid", "maps_to_dbsnp", "new_rsid", "merged_into_new_rsid",
  "bsgenome_objects", "has_rsid", "head", "merge_history", "only_b37", "history",
  "exists_only_on_grch37", "reason", "filter_callback", "no_dbsnp_entry", "logfile"))

snp_cols <- c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele", "rowid")
info_cols <- c("INFO", "N", "CaseN", "ControlN", "EAF")
stats_cols <- c("B", "Z", "OR", "P", "SE", info_cols)
valid_column_names <- c(snp_cols, stats_cols, info_cols)




#' Clean GWAS summary statistics
#'
#' @param tbl a tibble or filepath
#' @param bsgenome_objects use get_bsgenome()
#' @param outdir where to copy results after finished
#' @param logfile Write messages to a logfile?
#' @param name name of the output directory
#' @param ... arguments to other functions, such as study_n, build
#'
#' @return a tibble or NULL, depending on outdir
#' @export
#'
#' @examples \dontrun{
#' tidyGWAS(tbl = "my_dataframe", logfile = "true", name = "test_run", outdir = "gwas_sumstat_dir")
#' }
tidyGWAS <- function(tbl, bsgenome_objects, logfile=FALSE, name, outdir, keep_indels = TRUE, ...) {
  stopifnot("tbl is a mandatory argument" = !missing(tbl))

  # setup  ------------------------------------------------------------------
  if(missing(name)) name <- stringr::str_replace_all(Sys.time(), pattern = "[:-]", replacement = "_") |> stringr::str_replace(" ", "_")
  filepaths <- setup_pipeline_paths(name = name)
  if(logfile) {
    cli::cli_alert_info("Output is redirected to logfile: {.file {filepaths$logfile}}")
    withr::local_message_sink(filepaths$logfile)
    withr::local_output_sink(filepaths$logfile)
  }


  # start pipeline ----------------------------------------------------------
  cli::cli_h1("Running {.pkg tidyGWAS {packageVersion('tidyGWAS')}}")
  cli::cli_inform("Starting at {Sys.time()}")
  rows_start <- nrow(tbl)

  # write raw file
  cli::cli_inform("Writing out raw sumstats")
  if(!"rowid" %in% colnames(tbl)) tbl <- dplyr::mutate(tbl, rowid = as.integer(row.names(tbl)))
  data.table::fwrite(tbl, paste0(filepaths$base , "/raw_sumstats.gz"))

  # run initial checks
  struct <- initiate_struct(tbl = tbl, ..., filepaths = filepaths)

  # validate SNP identifiers  ----------------------------------------------
  struct <- validate_snps(struct, .filter_callback = make_callback(struct$filepaths$validate_snps))

  # Validate the stats columns ----------------------------------------------
  struct <- validate_stats(struct, .filter_callback = make_callback(struct$filepaths$validate_stats))


  # Checks against dbSNP ----------------------------------------------------
  if(!missing(bsgenome_objects)) struct$sumstat <- validate_with_dbsnp(struct, bsgenome_objects = bsgenome_objects, .filter_callback = make_callback(struct$filepaths$validate_with_dbsnp))

  if(keep_indels & nrow(struct$indels) > 0) {
    indel_struct <- vector("list")
    indel_struct$sumstat <- dplyr::select(struct$indels, dplyr::any_of(c("rowid", snp_cols)))
    indel_struct$stats <-   dplyr::select(struct$indels, dplyr::any_of(c("rowid", stats_cols, info_cols)))
    indel_struct$has_rsid <- FALSE
    indel_struct <- validate_snps(indel_struct, validate_alleles = FALSE, .filter_callback = make_callback(struct$filepaths$validate_snps_indels))
    indel_struct <- validate_stats(indel_struct, .filter_callback = make_callback(struct$filepaths$validate_stats_indels))

    # merge back into main struct
    cli::cli_alert_success("Indels have been merged back into main pipeline")
    struct$sumstat <- dplyr::bind_rows(struct$sumstat, indel_struct$sumstat)
    struct$stats <- dplyr::bind_rows(struct$stats, indel_struct$stats)
  }

  # main checks done --------------------------------------------------------
  main <- dplyr::inner_join(struct$sumstat, struct$stats, by = "rowid")


  # print info about end results ----------------------------------------------------------------
  cli::cli_h1("Finished tidyGWAS")
  cli::cli_alert_info("A total of {rows_start - nrow(main)} rows were removed. Started with: {rows_start} rows, ended with: {nrow(main)} rows")


  identify_removed_rows(dplyr::select(main,rowid), struct$filepaths)

  if(!missing(outdir)) {
    cleaned <- dplyr::select(main, dplyr::any_of(c("rowid", "CHR", "POS", "RSID","EffectAllele", "OtherAllele")), dplyr::everything())
    data.table::fwrite(cleaned, struct$filepaths$cleaned)
    file.copy("tidyGWAS_log.txt", paste0(struct$filepaths$base,"/tidyGWAS_log.txt"))
    file.copy(struct$filepaths$base, outdir, recursive = TRUE)
    return(NULL)
  } else {

    main
  }


}


validate_with_dbsnp <- function(struct, bsgenome_objects, .filter_callback, ...) {
  create_messages("validate_with_dbsnp")


  # existence of chr:pos or rsid decides which columns to repair
  # also checks if build has bee provided
  if(!struct$has_chr_pos & struct$has_rsid) {
    main_df <- repair_chr_pos(sumstat = struct$sumstat, bsgenome_objects = bsgenome_objects)

  } else if(struct$has_chr_pos & !struct$has_rsid) {
    # build can be passed with ...
    if(!is.null(struct$build)) {
      main_df <- repair_rsid(sumstat = struct$sumstat, bsgenome_objects = bsgenome_objects, build = struct$build)
    } else {
      main_df <- repair_rsid(sumstat = struct$sumstat, bsgenome_objects = bsgenome_objects)
    }

  } else if(struct$has_chr_pos & struct$has_rsid) {
    if(!is.null(struct$build)) {
      main_df <- verify_chr_pos_rsid(sumstat = struct$sumstat, bsgenome_objects = bsgenome_objects, build = struct$build)
    } else {
      main_df <- verify_chr_pos_rsid(sumstat = struct$sumstat, bsgenome_objects = bsgenome_objects)
    }
  }


  # check if there is a subset of SNPs detected as REF:ALT:CHR:POS  ---------
  if(!is.null(struct$without_rsid)) {
    cli::cli_h1("Found rows from RSID column which could be parsed as CHR:POS:REF:ALT or CHR:POS")
    cli::cli_inform("Validating the parsed data, and adding RSID")
    filepaths <- setup_pipeline_paths(rsid_subset = paste0(struct$filepaths$base, "/rsid_subset"))
    tmp <- initiate_struct(struct$without_rsid, filepaths = filepaths)
    tmp <- validate_snps(tmp, .filter_callback = make_callback(struct$filepaths$validate_snps))

    tmp <- repair_rsid(tmp$sumstat, bsgenome_objects = bsgenome_objects)
    main_df <- dplyr::bind_rows(main_df, tmp)
  }



  cli::cli_h2("Finished validation against dbSNP v.155")
  if(!missing(.filter_callback)) main_df <- .filter_callback(main_df)
  main_df


}




initiate_struct <- function(tbl, build, rs_merge_arch, filepaths, study_n,  ...) {

  # starting print
  cli::cli_h2("Performing initial checks on input data: ")
  cli::cli_ul()
  ulid <- cli::cli_ol()
  if(!missing(build)) cli::cli_alert_info("Using GRCh{build} as genome build")
  if(!missing(study_n)) cli::cli_alert_info("Using N={study_n} to impute N if N column is not present")
  if(missing(rs_merge_arch)) rs_merge_arch <- get_ref_data()


  # check input columns
  cli::cli_alert_info("Columns: {colnames(tbl)[colnames(tbl) %in% valid_column_names]}")
  if(length(colnames(tbl)[!colnames(tbl) %in% valid_column_names] > 0)) {
    cli::cli_alert_danger(c("{.strong Removed columns: }",
      "{colnames(tbl)[!colnames(tbl) %in% valid_column_names]}"))
  }
  # check that either RSID or CHR:POS is present
  stopifnot("Requires either RSID or CHR:POS" = c("RSID") %in% colnames(tbl) | c("CHR", "POS") %in% colnames(tbl))
  stopifnot("Requires EffectAllele and OtherAllele" = all(c("EffectAllele", "OtherAllele") %in% colnames(tbl)))

  # setup filepaths that will be used, and remove unwanted columns
  struct <- vector("list")
  struct$filepaths <- filepaths
  tbl <- dplyr::select(tbl,dplyr::any_of(valid_column_names))

  # add rowid if it doesnt exist
  if(!"rowid" %in% colnames(tbl)) tbl <- dplyr::mutate(tbl, rowid = as.integer(row.names(tbl)))

  # add build if it exists
  if(missing(build)) struct[["build"]] <- NULL else struct[["build"]] <- build

  # add N if CaseN and ControlN exists
  if(all(c("CaseN", "ControlN") %in% colnames(tbl))) tbl$N <- (tbl$CaseN + tbl$ControlN)

  # impute_n if argument supplied and infor if N is missing
  if(!missing(study_n)) tbl <- dplyr::mutate(tbl, N = {{ study_n }})
  if(!"N" %in% colnames(tbl)) cli::cli_alert_danger("Found no N column, and no study_n was supplied. It is highly recommended to supply a value for N, as many downstream GWAS applications rely on this information")

  # add flags for which identifiers exist
  struct[["has_chr_pos"]] <- ifelse("CHR"  %in% colnames(tbl) & "POS" %in% colnames(tbl), TRUE, FALSE)
  struct[["has_rsid"]] <- ifelse("RSID" %in% colnames(tbl), TRUE, FALSE)


  # remove rows with NA -----------------------------------------------------

  tmp <- tidyr::drop_na(tbl, -dplyr::any_of(c("CHR", "POS", "RSID")))
  struct$na_rows <- dplyr::anti_join(tbl, tmp, by = "rowid")
  create_messages(func = "drop_na", struct = struct)

  # handle duplicates -------------------------------------------------------

  cli::cli_h3("Looking for duplicate rows")
  cli::cli_inform("If possible, selects the row with smallest pvalue in each duplication unit")
  if("P" %in% colnames(tmp)) tmp <- dplyr::arrange(tmp, .data[["P"]])

  if(struct$has_chr_pos) {
    id <- "CHR_POS_REF_ALT"
    cols_to_use <- c("CHR", "POS", "EffectAllele", "OtherAllele")
  } else {
    id <- "RSID_REF_ALT"
    cols_to_use <- c("RSID", "EffectAllele", "OtherAllele")

  }

  cli::cli_alert_info("Using {id} as id")
  # use distinct to remove duplications
  no_dups <- dplyr::distinct(tmp, dplyr::pick(dplyr::all_of(cols_to_use)), .keep_all = TRUE)
  removed <- dplyr::anti_join(tmp, no_dups, by = "rowid")
  tmp <- no_dups
  if(nrow(removed) > 0) cli::cli_alert_info("Removed {nrow(removed)} rows flagged as duplications")



  # remove indels from main pipeline ------------------------------------------
  create_messages("indels")
  tmp <- flag_indels(tmp)
  struct$indels <-  dplyr::select(dplyr::filter(tmp,  .data[["indel"]]), -indel)
  tmp <-     dplyr::select(dplyr::filter(tmp, !.data[["indel"]]), -indel)


  # if possible, update RSID ------------------------------------------------

  if("RSID" %in% colnames(tmp)) {
    cli::cli_h3("Updating RSIDs that have been merged using RsMergeArch")
    rsid_info <- flag_rsid_history(tmp, rs_merge_arch = rs_merge_arch)
    tmp <- dplyr::inner_join(dplyr::select(tbl, -RSID), rsid_info, by = "rowid")
    cli::cli_alert_success("Updated RSID in {sum(!is.na(rsid_info$old_RSID))} rows")

  } else {
    rsid_info <- dplyr::tibble(old_RSID = "X", .rows = 0)
  }


  # summarise what has been doen --------------------------------------------
  cli::cli_h2("Finished initial checks")
  cli::cli_ul()


  if(sum(!is.na(rsid_info$old_RSID)) > 0) {
    cli::cli_li("{sum(!is.na(rsid_info$old_RSID))} rows with updated RSID: {.file {struct$filepaths$updated_rsid}}")
    data.table::fwrite(dplyr::filter(rsid_info, !is.na(old_RSID)), struct$filepaths$updated_rsid, sep = "\t")
  }

  if(nrow(removed) > 0) {
    cli::cli_li("{nrow(removed)} rows removed because duplicates: {.file {struct$filepaths$duplicates}}")
    data.table::fwrite(removed, struct$filepaths$duplicates)
  }

  # split into statistics and snp identifiers -------------------------------

  struct$sumstat <- dplyr::select(tmp, dplyr::any_of(c("rowid", snp_cols)))
  struct$stats <-   dplyr::select(tmp, dplyr::any_of(c("rowid", stats_cols, info_cols)))


  # return results ----------------------------------------------------------
  struct

}

validate_stats <- function(struct, ..., .filter_callback) {
  create_messages("validate_stats", struct$stats)
  tbl <- struct$stats
  # check for P = 0 and imputation of Z
  warning_messages_stats(tbl)




  # run validators ----------------------------------------------------------

  if("EAF" %in% colnames(tbl)) tbl <- validate_eaf(tbl)
  if("SE"  %in% colnames(tbl)) tbl <- validate_se(tbl)
  if("P"   %in% colnames(tbl)) tbl <- validate_P(tbl, ...)
  if("B"   %in% colnames(tbl)) tbl <- validate_b(tbl)
  if("Z"   %in% colnames(tbl)) tbl <- validate_z(tbl)
  if("N"   %in% colnames(tbl)) tbl <- validate_n(tbl)


  if(!missing(.filter_callback)) tbl <- .filter_callback(tbl)


  struct$stats <- tbl
  struct
}



validate_snps <- function(struct, .filter_callback, validate_alleles=TRUE) {

  create_messages("validate_snps")

  if(struct$has_rsid) {

    # unpack validation of RSIDS
    # need to process any CHR:POS:REF:ALT separately
    validated_rsid <- validate_rsid(struct$sumstat)
    struct$sumstat <- validated_rsid$data
    struct$without_rsid <- validated_rsid$chr_pos
    struct$failed <-  validated_rsid$failed
  }

  if(!is.null(struct$failed)) {
    if(nrow(struct$failed > 0)) data.table::fwrite(struct$failed, struct$filepaths$failed_rsid_parse)
  }



  if("CHR" %in% colnames(struct$sumstat)) struct$sumstat <- validate_chr(struct$sumstat)
  if("POS" %in% colnames(struct$sumstat)) struct$sumstat <- validate_pos(struct$sumstat)
  if(validate_alleles) struct$sumstat <- validate_ea_oa(struct$sumstat)

  # use callback
  if(!missing(.filter_callback)) struct$sumstat <- .filter_callback(struct$sumstat)


  struct

}


warning_messages_stats <- function(tbl) {
  n_p_0 <- dplyr::filter(tbl, P == 0) |> nrow()
  z_miss <- all(c("B", "P") %in% colnames(tbl)) & !"Z" %in% colnames(tbl)
  if(n_p_0 & z_miss) cli::cli_alert_danger("WARNING: Found {n_p_0} rows with P = 0 and missing Z score. These will not be correctly imputed")
}


setup_pipeline_paths <- function(name, rsid_subset) {

  if(missing(rsid_subset)) {
    if(!dir.exists(tempdir())) dir.create(tempdir())
    workdir <- paste(tempdir(), name, name, sep = "/")

  } else {
    workdir <- rsid_subset
  }

  pipeline_info <- paste(workdir,"pipeline_info", sep = "/")
  indel_dir <- paste(pipeline_info,"indels", sep = "/")
  suppressWarnings(dir.create(pipeline_info, recursive = TRUE))
  if(!dir.exists(indel_dir)) dir.create(indel_dir)



  # name/name/:
  # /pipeline_info/
  # /raw
  # /cleaned
  #


  start_ids <-           paste(pipeline_info, "start_ids.tsv.gz", sep = "/")
  duplicates <-          paste(pipeline_info, "duplicates.log.gz", sep = "/")
  validate_stats <-      paste(pipeline_info, "validate_stats.log.gz", sep = "/")
  validate_snps <-       paste(pipeline_info, "validate_snps.log.gz", sep = "/")
  validate_with_dbsnp <- paste(pipeline_info, "validate_with_dbsnp.log.gz", sep = "/")
  updated_rsid <-        paste(pipeline_info, "updated_rsid.log.gz", sep = "/")
  rows_with_na <-        paste(pipeline_info, "rows_with_na.log.gz", sep = "/")
  failed_rsid_parse <-   paste(pipeline_info, "failed_rsid_parse.log.gz", sep = "/")


  indels <-              paste(indel_dir, "indels.log.gz", sep = "/")
  validate_stats_indels <-      paste(indel_dir, "validate_stats.log.gz", sep = "/")
  validate_snps_indels <-       paste(indel_dir, "validate_snps.log.gz", sep = "/")

  invalid_rsid <-  paste(pipeline_info, "invalid_rsid", sep = "/")
  if(!dir.exists(invalid_rsid)) dir.create(invalid_rsid)
  invalid_rsid <-  paste0(workdir, "/pipeline_info/invalid_rsid")




  list(
    "base" = workdir,
    "logfile" = paste0(workdir, "/tidyGWAS_logfile.txt"),
    "cleaned" = paste(workdir, "cleaned_GRCh38.gz", sep = "/"),
    "start_ids" = start_ids,
    "duplicates" = duplicates,
    "validate_stats" = validate_stats,
    "validate_snps" = validate_snps,
    "validate_with_dbsnp" = validate_with_dbsnp,
    "updated_rsid"= updated_rsid,
    "failed_rsid_parse" = failed_rsid_parse,
    "indels" =                     paste(indel_dir, "indels.log.gz", sep = "/"),
    "validate_stats_indels" =      paste(indel_dir, "validate_stats.log.gz", sep = "/"),
    "validate_snps_indels" =       paste(indel_dir, "validate_snps.log.gz", sep = "/"),
    "invalid_rsid" = invalid_rsid,
    "rows_with_na" = rows_with_na
  )


}

create_messages <- function(func,tbl, struct) {
  if(func== "validate_snps") {
    cli::cli_h2("Starting validation of CHR,POS, RSID, EffectAllele and OtherAllele")

  # validate_with_dbnsp -----------------------------------------------------


  } else if(func == "validate_with_dbsnp") {
    cli::cli_h2("Validating sumstats using dbSNP")
    cli::cli_ol()
    cli::cli_li("Repair missing CHR, POS or RSID")
    cli::cli_li("Remove rows where REF/ALT in dbSNP is not compatible with EffectAllele / OtherAllele")
    cli::cli_li("Remove rows where rsID does not match any entry in dbSNP v.155")

  # validate_stats ----------------------------------------------------------


  } else if(func == "validate_stats") {
    cols_in_tbl <- colnames(tbl)[colnames(tbl) %in% stats_cols]

    cli::cli_h2("Validating statistics columns: {cols_in_tbl}")
  }

  # drop MISSING ----------------------------------------------------------------


  else if(func=="drop_na") {
    cli::cli_h3("Examining data for missing values across all columns")
    cli::cli_ul("Using tidyr::drop_na()")
    if(nrow(struct$na_rows > 0)){
      cli::cli_alert("Found {nrow(struct$na_rows)} rows with missing values. These are removed. Printing the first 5 rows")
      cli::cat_print(dplyr::slice(struct$na_rows, 1:5))
      data.table::fwrite(struct$na_rows, struct$filepaths$rows_with_na, sep = "\t")
    } else {
      cli::cli_alert_success("Found no rows with missing values")
    }


  # INDELS ------------------------------------------------------------------


  } else if(func == "indels") {
    cli::cli_h3("Checking for insertions/deletions ('indels') using: ")
    cli::cli_ol(c(
      "EffectAllele or OtherAllele, character length > 1: A vs AA",
      "EffectAllele or OtherAllele coded as 'D', 'I', or 'R'"
    ))



  }


}


identify_removed_rows <- function(finished, filepaths) {
  start <- data.table::fread(paste0(filepaths$base , "/raw_sumstats.gz"), select = "rowid")

  removed_rows <- dplyr::anti_join(start, finished, by = "rowid")

  # remove the file with all rowids
  files_in_dir <- list.files(paste0(filepaths$base, "/pipeline_info/"), full.names = TRUE, pattern = "*.gz")
  files_in_dir <- files_in_dir[!stringr::str_detect(files_in_dir, "start_ids.tsv.gz")]
  files_in_dir <- files_in_dir[!stringr::str_detect(files_in_dir, "updated_rsid.log.gz")]


  removed_rows_with_flags <-
    files_in_dir |>
    purrr::set_names(base::basename) |>
    purrr::map(data.table::fread) |>
    purrr::map(\(x) dplyr::select(x, rowid)) |>
    purrr::list_rbind(names_to = "reason") |>
    dplyr::select(rowid, reason) |>
    dplyr::mutate(reason = stringr::str_remove(reason, ".log.gz"))

  if(sum(!removed_rows$rowid %in% removed_rows_with_flags$rowid)) {
    cli::cli_alert_danger("WARNING: Could not track why some rows were removed. This is very surprising")
  }

  breakdown <-
    removed_rows_with_flags |>
    dplyr::semi_join(removed_rows, by = "rowid") |>
    dplyr::count(reason)

  vec <- c(breakdown$n)
  names(vec) <- breakdown$reason
  cli::cli_h3("Listing final breakdown of removed rows: ")
  cli::cli_dl(vec)




}

# Suppress R CMD check note
#' @importFrom R.utils as.character.binmode
NULL
