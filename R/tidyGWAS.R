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
#' @param use_dbsnp use dbSNP to apply filters? Default is TRUE
#' @param outdir Where should results be saved after a succesful run? Default is tempdir()
#' @param logfile Direct messages to a logfile? Default is FALSE
#' @param name name of the output directory. Default is a concotonated call to Sys.time()
#' @param keep_indels Should indels be kept? Default is TRUE
#' @param verbose Explain filters in detail? Default is FALSE.
#' @param log_on_err Optional. Can pass a filepath to copy the logfile to when the function exists.
#' This can be very useful if running not interactively, and want to make sure
#' the log file exists even if the function errors.
#' @param ... arguments that will be passed to other functions: Currently supports
#' study_n, build, bsgenome_objects, rs_merge_arch
#'
#' @return a tibble or NULL, depending on outdir
#' @export
#'
#' @examples \dontrun{
#' tidyGWAS(tbl = "my_dataframe", logfile = "true", name = "test_run", outdir = "gwas_sumstat_dir")
#' }
tidyGWAS <- function(
    tbl,
    use_dbsnp = TRUE,
    name = stringr::str_replace_all(date(), pattern = c(" "="_", ":"="_")),
    outdir,
    bsgenome_objects,
    logfile=FALSE,
    log_on_err="tidyGWAS_log.txt",
    keep_indels = TRUE,
    verbose = FALSE,
    ...
    ) {

  # parse data and input checks ------------------------------------------------

  stopifnot("tbl is a mandatory argument" = !missing(tbl))
  stopifnot("tbl is not character or data.frame" = any(c("data.frame", "character") %in% class(tbl)))
  filepaths <- setup_pipeline_paths(name = name)

  if("character" %in% class(tbl)) {
    tbl <- data.table::fread(tbl)
  } else if("data.frame" %in% class(tbl)) {
    tbl <- dplyr::tibble(tbl)
  }
  rows_start <- nrow(tbl)
  start_time <- Sys.time()

  # use withr to create a sink that cancels on exit
  if(logfile) {
    cli::cli_alert_info("Output is redirected to logfile: {.file {filepaths$logfile}}")
    withr::local_message_sink(filepaths$logfile)
    withr::local_output_sink(filepaths$logfile)
    if(!missing(log_on_err)) on.exit(file.copy(filepaths$logfile, log_on_err), add=TRUE)
  }

  # check for bsgenome_object
  if(missing(bsgenome_objects)) bsgenome_objects <- get_bsgenome()


  # welcome message ----------------------------------------------------------

  cli::cli_h1("Running {.pkg tidyGWAS {packageVersion('tidyGWAS')}}")
  cli::cli_inform("Starting at {start_time}")
  cli::cli_alert_info("Saving all files during execution to {.file {filepaths$base}}")
  if(!missing(outdir)) cli::cli_alert_info("Files will be copied to {.file {paste(outdir,name, sep = '/')}} when finished")



  # start the pipeline ------------------------------------------------------

  # run initial checks
  struct <- initiate_struct(tbl = tbl, filepaths = filepaths, verbose = verbose, ...)

  # validate sumstats
  struct <- validate_sumstat(struct, verbose = verbose)

  # Validate with dbSNP
  if(use_dbsnp) {
    struct$sumstat <- validate_with_dbsnp(
      struct,
      bsgenome_objects = bsgenome_objects
      )
  }




  # merge indels back?
  if(keep_indels & nrow(struct$indels) > 0) {
    indel_struct <- vector("list")
    indel_struct$sumstat <- struct$indels
    indel_struct$filepaths <- struct$filepaths

    # We can't validate RSIDs for indels, as indels are not in the BSgenome dbsnp version. Might be an upcoming feature.
    indel_struct$has_rsid <- FALSE
    indel_struct <- validate_sumstat(indel_struct, remove_cols = c("EffectAllele","OtherAllele"), verbose = verbose)

    # merge back into main struct
    cli::cli_alert_success("Indels have been merged back into main pipeline")
    struct$sumstat <- dplyr::bind_rows(struct$sumstat, indel_struct$sumstat)
  }
  main <- struct$sumstat

  # Finished! wrap up ----------------------------------------------

  cli::cli_h1("Finished tidyGWAS")
  cli::cli_alert_info("A total of {rows_start - nrow(main)} rows were removed")
  identify_removed_rows(dplyr::select(main,rowid), struct$filepaths)




  # return cleaned or write out? ----------------------------------------------


  if(!missing(outdir)) {
    cleaned <- dplyr::select(main, dplyr::any_of(c("rowid", "CHR", "POS", "RSID","EffectAllele", "OtherAllele")), dplyr::everything())
    data.table::fwrite(cleaned, struct$filepaths$cleaned)
    file.copy("tidyGWAS_log.txt", paste0(struct$filepaths$base,"/tidyGWAS_log.txt"))
    file.copy(struct$filepaths$base, outdir, recursive = TRUE)
    cli::cli_li("Total running time: {prettyunits::pretty_dt(end_time - start_time)}")


  } else {
    end_time <- Sys.time()
    cli::cli_li("Total running time: {prettyunits::pretty_dt(end_time - start_time)}")
    main
  }



}


# -------------------------------------------------------------------------


validate_with_dbsnp <- function(struct, bsgenome_objects) {
  create_messages("validate_with_dbsnp")
  .filter_callback = make_callback(struct$filepaths$validate_with_dbsnp)

  # existence of chr:pos or rsid decides which columns to repair
  if(!struct$has_chr_pos & struct$has_rsid) {
    main_df <- repair_chr_pos(sumstat = struct$sumstat, bsgenome_objects = bsgenome_objects)

  } else if(struct$has_chr_pos & !struct$has_rsid) {
    # nested if else to pass build if it was given in tidyGWAS
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


  # check if there is a subset of SNPs with missing rsid that needs to be repaired
  if(!is.null(struct$without_rsid)) {
    if(nrow(struct$without_rsid) > 0) {
      tmp <- repair_rsid(struct$without_rsid, bsgenome_objects = bsgenome_objects)
      main_df <- dplyr::bind_rows(main_df, tmp)
    }
  }



  cli::cli_h2("Finished validation against dbSNP v.155")
  if(!missing(.filter_callback)) main_df <- .filter_callback(main_df)
  main_df


}


# -------------------------------------------------------------------------




initiate_struct <- function(tbl, filepaths, verbose=FALSE, rs_merge_arch, study_n, ...) {


  # -------------------------------------------------------------------------



  cli::cli_h2("Performing initial checks on input data: ")
  if(!missing(study_n)) cli::cli_alert_info("Using N={study_n} to impute N if N column is not present")
  if(missing(rs_merge_arch)) rs_merge_arch <- get_ref_data()


  # check input columns
  cli::cli_alert_info("Columns: {colnames(tbl)[colnames(tbl) %in% valid_column_names]}")
  if(length(colnames(tbl)[!colnames(tbl) %in% valid_column_names] > 0)) {
    cli::cli_alert_danger("{.strong Removed columns:  {colnames(tbl)[!colnames(tbl) %in% valid_column_names]}}")
  }
  # check that either RSID or CHR:POS is present
  stopifnot("Requires either RSID or CHR:POS" = c("RSID") %in% colnames(tbl) | c("CHR", "POS") %in% colnames(tbl))
  stopifnot("Requires EffectAllele and OtherAllele" = all(c("EffectAllele", "OtherAllele") %in% colnames(tbl)))

  # write raw file
  tbl <- dplyr::tibble(tbl)
  cli::cli_inform("Keeping track of rows by writing out a rowindex file of raw input file")
  if(!"rowid" %in% colnames(tbl)) tbl <- dplyr::mutate(tbl, rowid = as.integer(row.names(tbl)))
  data.table::fwrite(dplyr::select(tbl, rowid), paste0(filepaths$base , "/raw_sumstats.gz"))

  # setup filepaths that will be used, and remove unwanted columns
  struct <- vector("list")
  struct$filepaths <- filepaths
  tbl <- dplyr::select(tbl,dplyr::any_of(valid_column_names))



  # add N if CaseN and ControlN exists
  if(all(c("CaseN", "ControlN") %in% colnames(tbl))) tbl$N <- (tbl$CaseN + tbl$ControlN)

  # impute_n if argument supplied and infor if N is missing
  if(!missing(study_n)) tbl <- dplyr::mutate(tbl, N = {{ study_n }})
  if(!"N" %in% colnames(tbl)) cli::cli_alert_danger("Found no N column, and no study_n was supplied. It is highly recommended to supply a value for N, as many downstream GWAS applications rely on this information")

  # add flags for which identifiers exist
  struct[["has_chr_pos"]] <- ifelse("CHR"  %in% colnames(tbl) & "POS" %in% colnames(tbl), TRUE, FALSE)
  struct[["has_rsid"]] <- ifelse("RSID" %in% colnames(tbl), TRUE, FALSE)


  # 1) First filter - remove rows with NA ---------------------------------------

  cli::cli_ul()
  cli::cli_li("Scanning for NAs with tidyr::drop_na()")

  tmp <- tidyr::drop_na(tbl, -dplyr::any_of(c("CHR", "POS", "RSID")))
  struct$na_rows <- dplyr::anti_join(tbl, tmp, by = "rowid")

  if(nrow(struct$na_rows) > 0) {
    cli::cli_li("Found {nrow(struct$na_rows)} rows with missing values. These are removed. Printing the first 5 rows")
    cli::cat_print(dplyr::slice(struct$na_rows, 1:5), file = stderr())
    data.table::fwrite(struct$na_rows, struct$filepaths$rows_with_na, sep = "\t")
  } else {
    cli::cli_alert_success("Found no rows with missing values")
  }


  # 2) if possible, update RSID ------------------------------------------------

  if(struct$has_rsid) {

    rsid_info <- flag_rsid_history(tmp, rs_merge_arch = rs_merge_arch)
    tmp <- dplyr::inner_join(dplyr::select(tbl, -RSID), rsid_info, by = "rowid") |>
      dplyr::select(-old_RSID)

    cli::cli_h3("Updating RSIDs that have been merged using RsMergeArch")
    cli::cli_li("{sum(!is.na(rsid_info$old_RSID))} rows with updated RSID: {.file {struct$filepaths$updated_rsid}}")
    data.table::fwrite(dplyr::filter(rsid_info, !is.na(old_RSID)), struct$filepaths$updated_rsid, sep = "\t")

  } else {
    rsid_info <- dplyr::tibble(old_RSID = "X", .rows = 0)
  }



  # 3) handle duplicates -------------------------------------------------------


  if("P" %in% colnames(tmp)) tmp <- dplyr::arrange(tmp, .data[["P"]])

  if(struct$has_chr_pos) {
    id <- "CHR_POS_REF_ALT"
    cols_to_use <- c("CHR", "POS", "EffectAllele", "OtherAllele")
  } else {
    id <- "RSID_REF_ALT"
    cols_to_use <- c("RSID", "EffectAllele", "OtherAllele")
  }

  cli::cli_li("Scanning for duplicates")
  cli::cli_li("Using {id} as id")

  # use distinct to remove duplications
  no_dups <- dplyr::distinct(tmp, dplyr::pick(dplyr::all_of(cols_to_use)), .keep_all = TRUE)
  removed <- dplyr::anti_join(tmp, no_dups, by = "rowid")
  tmp <- no_dups

  if(nrow(removed) > 0) {
    cli::cli_alert_danger("Removed {nrow(removed)} rows flagged as duplications, based on {cols_to_use}")
    cli::cli_li("{nrow(removed)} rows removed because duplicates: {.file {struct$filepaths$duplicates}}")
    data.table::fwrite(removed, struct$filepaths$duplicates)
  }


  # 4) remove indels from main pipeline ------------------------------------------

  if(verbose) create_messages("indels")
  tmp <- flag_indels(tmp)
  struct$indels <-  dplyr::select(dplyr::filter(tmp,  .data[["indel"]]), -indel)
  tmp <-            dplyr::select(dplyr::filter(tmp, !.data[["indel"]]), -indel)



  # return results ----------------------------------------------------------
  cli::cli_h2("Finished initial checks")
  struct$sumstat <- tmp
  struct

}



# -------------------------------------------------------------------------


validate_sumstat <- function(struct, remove_cols="", verbose=FALSE, filter_rows=TRUE) {

  stopifnot("remove_cols can only be a character vector" = is.character(remove_cols))
  # columns which has a validtor
  impl_validators <- c("CHR", "POS", "EffectAllele", "OtherAllele","EAF", "SE", "P", "B", "Z", "N")
  # columns in remove_cols should not be validated
  impl_validators <- impl_validators[!impl_validators %in% remove_cols]

  # RSID validation ---------------------------------------------------------
  # a bit more complex than other columns
  # need to handle rows where RSID is actually CHR:POS:REF:ALT or similar


  if(struct$has_rsid) {
    validated_rsid <- validate_rsid(struct$sumstat, verbose = verbose)
    # after calling validation, rows with valid RSIDs are kept in main pipeline
    struct$sumstat <- validated_rsid$data |> dplyr::filter(!invalid_rsid) |> dplyr::select(-invalid_rsid)
    # succesful parses follow another path, and will have their rsid repair in validate_with_dbsnp
    struct$without_rsid <- validated_rsid$chr_pos
    # failed parsed are removed
    struct$failed <-  validated_rsid$failed
    if(!is.null(struct$failed)) if(nrow(struct$failed > 0)) data.table::fwrite(struct$failed, struct$filepaths$failed_rsid_parse)
    if(is.null(struct$without_rsid)) struct$without_rsid <- dplyr::tibble(.rows = 0)
    # validate the subset without RSID (CHR,POS, EA, OA)
    if(nrow(struct$without_rsid) != 0) {
      cli::cli_h3("Running validation of rows without a valid RSID")
      cols_to_run <- colnames(struct$without_rsid)[colnames(struct$without_rsid) %in% impl_validators]
      for(c in cols_to_run) struct$without_rsid <- validate_columns(struct$without_rsid, col = c, verbose = FALSE)
      tbl <- repair_stats(tbl)


      cb_new <- make_callback(struct$filepaths$validate_snps)
      struct$without_rsid <- cb_new(struct$without_rsid)
    }
  }




  # run validators on main data------------------------------------------------

  tbl <- struct$sumstat
  warning_messages_stats(tbl)

  cli::cli_h3("Running validation for other rows")
  cols_in_sumstat <- colnames(tbl)[colnames(tbl) %in% impl_validators]
  for(c in cols_in_sumstat) tbl <- validate_columns(tbl, col = c, verbose = verbose)

  tbl <- repair_stats(tbl)
  # write failed rows to disk if passed
  if(filter_rows) {
    # if some rows have already been removed in validate SNPs, need to make sure to not overwrite file
    filtering_function <- make_callback(outpath = struct$filepaths$validate_snps)
    tbl <- filtering_function(tbl)
  }
  struct$sumstat <- tbl


  # finished ----------------------------------------------------------------


  struct

}


# -------------------------------------------------------------------------


warning_messages_stats <- function(tbl) {
  n_p_0 <- dplyr::filter(tbl, P == 0) |> nrow()
  z_miss <- all(c("B", "P") %in% colnames(tbl)) & !"Z" %in% colnames(tbl)
  if(n_p_0 & z_miss) cli::cli_alert_danger("WARNING: Found {n_p_0} rows with P = 0 and missing Z score. These will not be correctly imputed")
}


# -------------------------------------------------------------------------


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

# -------------------------------------------------------------------------


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


# -------------------------------------------------------------------------


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

# exit <- function(main, outdir) {
#   if(!missing(outdir)) {
#     cleaned <- dplyr::select(main, dplyr::any_of(c("rowid", "CHR", "POS", "RSID","EffectAllele", "OtherAllele")), dplyr::everything())
#     data.table::fwrite(cleaned, struct$filepaths$cleaned)
#     file.copy("tidyGWAS_log.txt", paste0(struct$filepaths$base,"/tidyGWAS_log.txt"))
#     file.copy(struct$filepaths$base, outdir, recursive = TRUE)
#     return(NULL)
#   } else {
#
#     main
#   }
# }
