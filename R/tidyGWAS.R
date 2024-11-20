utils::globalVariables(c(
  "chr_pos", "n", "n_chr_pos", "dup_rsid",
  ".data", "invalid_rsid",  "new_rsid",
  "has_rsid", "head", "CHR_37", "POS_37",
   "reason", "filter_callback", "no_dbsnp_entry", "logfile",
  "alt_allele"))

snp_cols <- c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele", "rowid")
info_cols <- c("INFO", "N", "CaseN", "ControlN","EffectiveN", "EAF")
stats_cols <- c("B", "Z", "OR", "P", "SE", info_cols)
valid_column_names <- c(snp_cols, stats_cols, info_cols)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------




#' Execute validation and quality control of GWAS summmary statistics
#'
#' @description
#' `tidyGWAS()` performs a set of validations on input colummns, repairs missing
#' columns, and can add missing CHR/POS or RSID. In addition, CHR and POS is
#' standardised to GRCh38, with coordinates on GRCh37 added in as well.
#'
#' Briefly, `tidyGWAS()` updates RSID if possible using the
#' [refsnp-merged](https://ftp.ncbi.nih.gov/snp/latest_release/JSON/) file
#' from dbSNP. Each inputed column is then validated and coerced to
#' the correct type.
#'
#' If statistis such as `P`, `B` are missing, `tidyGWAS()` will attempt to impute
#' them if possible using [repair_stats()]
#'
#' Standard column names are assumed, BEFORE inputting into the function. This is
#' a deliberate decision as automatic parsing of some important column names
#' can be ambiguous For example, in some sumstats, A1 referes to effect allele,
#' while other formats use A1 as non-effect allele.
#'
#' @param tbl a `data.frame` or `character()` vector
#' @param dbsnp_path filepath to the dbSNP155 directory (untarred dbSNP155.tar)
#' @param ... pass additional arguments to [arrow::read_delim_arrow()], if tbl is a filepath.
#' @param output_dir filepath to a folder where tidyGWAS output will be stored.
#'   The folder should not yet exist. Note that the default argument is `tempfile()`,
#'   meaning that tidyGWAS output will not be saved by default over R sessions.
#' @param output_format How should the finished cleaned file be saved?
#'  * 'csv' corresponds to [arrow::write_csv_arrow()]
#'  * 'parquet' corresponds to [arrow::write_parquet()]
#'  * 'hivestyle' corresponds to [arrow::write_dataset()] split by `CHR`
#'
#' @param logfile Should messages be redirected to a logfile?
#' @param column_names a named list of column names:
#' `list(RSID = "SNP", POS = "BP")`
#' @param CaseN manually input number of cases
#' @param ControlN manually input number of controls
#' @param N manually input sample size
#' @param impute_freq Should allele frequency be imputed if it's missing?
#'  Provide a filepath to a .parquet file with columns RSID, EffectAllele, OtherAllele, EAF.
#'  where EAF corresponds to frequency of the EffectAllele.
#' @param impute_n Should N be imputed if it's missing? see discussion in:
#'  https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics
#' @param allow_duplications Should duplicated variants be allowed? Useful if the munged sumstats are QTL sumstats
#' @param build If you are sure of what genome build ('37' or '38'), can be used to skip [infer_build()] and speed up computation
#' @param convert_p What value should be used for when P-value has been rounded to 0?
#' @param indel_strategy Should indels be kept or removed?
#' @param repair_cols Should any missing statistical columns be repaired if possible? calls [repair_stats()] if TRUE
#' @param default_build If only RSID exists, the build cannot be inferred. Nonetheless,
#' tidyGWAS applies a filter on incompatible alleles with GRCh37/38. In such a case,
#' tidyGWAS needs to decide on which reference genome to compare alleles with.
#'
#' @return a [dplyr::tibble()]
#'
#'
#' @export
#'
#' @examples \dontrun{
#' tidyGWAS(tbl = "my_dataframe", logfile = "true", name = "test_run", outdir = "gwas_sumstat_dir")
#' }
tidyGWAS <- function(
    tbl,
    dbsnp_path,
    ...,
    column_names,
    output_format = c("hivestyle","parquet", "csv"),
    output_dir = tempfile(),
    CaseN = NULL,
    ControlN = NULL,
    N = NULL,
    impute_freq = NULL,
    impute_n = FALSE,
    allow_duplications = FALSE,
    build = c("NA","37", "38"),
    default_build = c("37", "38"),
    indel_strategy = c("keep", "remove"),
    convert_p = 2.225074e-308,
    repair_cols = TRUE,
    logfile = FALSE
    ) {



  # parse arguments ---------------------------------------------------------
  rlang::check_required(tbl)
  rlang::check_required(dbsnp_path)
  if(!is.null(CaseN)) stopifnot(rlang::is_scalar_integerish(CaseN))
  if(!is.null(ControlN))  stopifnot(rlang::is_scalar_integerish(ControlN))
  if(!is.null(N)) stopifnot(rlang::is_scalar_integerish(N))
  if(!is.null(impute_freq)) stopifnot(rlang::is_scalar_character(impute_freq))
  if(!is.null(impute_freq)) stopifnot(file.exists(impute_freq))
  stopifnot(rlang::is_scalar_double(convert_p))
  stopifnot(rlang::is_bool(repair_cols))
  stopifnot(rlang::is_bool(logfile))
  stopifnot(rlang::is_bool(allow_duplications))
  stopifnot("logfile can only be TRUE or FALSE"= rlang::is_bool(logfile))
  stopifnot("The `output_dir` specified already exists" = !dir.exists(output_dir))
  stopifnot("the filepath for dbSNP does not exist" = file.exists(dbsnp_path))
  output_format <-  rlang::arg_match(output_format)
  build <-          rlang::arg_match(build)
  default_build <- rlang::arg_match(default_build)
  indel_strategy <- rlang::arg_match(indel_strategy)


  # welcome message ----------------------------------------------------------
  start_time <- Sys.time()
  cli::cli_h1("Running {.pkg tidyGWAS {packageVersion('tidyGWAS')}}")
  cli::cli_inform("Starting at {start_time}")


  # check input data.frame --------------------------------------------------

  tbl <- parse_tbl(tbl, ...)
  filename <- tbl$filename
  md5 <- tbl$md5
  # has to be last, overwriting tbl here
  tbl <- tbl$tbl

  rows_start <- nrow(tbl)
  filepaths <- setup_pipeline_paths(output_dir, filename)

  # setup logging -----------------------------------------------------------
  if(isTRUE(logfile)) {
    cli::cli_alert_info("Output is redirected to logfile: {.file {filepaths$logfile}}")
    withr::local_message_sink(filepaths$logfile)
    withr::local_output_sink(filepaths$logfile)
  }
  cli::cli_inform("with {rows_start} rows in input data.frame")
  cli::cli_alert_info("Saving output in folder: {.file {filepaths$base}}")

  # The sumstats are saved without any edits, to not loose information
  arrow::write_parquet(tbl,  filepaths$raw_sumstats)


  # start of pipeline ----------------------------------------------------------
  # 0) formatting --------------------------------------------------------------


  if(!missing(column_names)) tbl <- update_column_names(tbl, column_names, CaseN = CaseN, ControlN=ControlN, N=N)

  tbl <- select_correct_columns(tbl)


  # 1) detect indels ----------------------------------------------------------
  cli::cli_h2("1) Scanning for indels")

  tbl <- remove_rows_with_na(tbl, c("EffectAllele", "OtherAllele"), filepaths$removed_missing_alleles)
  tbl <- detect_indels(tbl, indel_strategy = indel_strategy, filepaths = filepaths, convert_p = convert_p)

  indels <- tbl$indels
  if(nrow(indels) == 0) indels <- NULL
  tbl <- tbl$main |> dplyr::mutate(indel = FALSE)


  # 2 Drop na -----------------------------------------------------------------
  if(all(c("RSID", "CHR", "POS") %in% colnames(tbl))) {
    tbl$RSID <- NULL
    columns <- c("CHR", "POS")
  } else if("RSID" %in% colnames(tbl) & !all(c("CHR", "POS") %in% colnames(tbl))) {
    columns <- c("RSID")
  } else {
    columns <- c("CHR", "POS")
  }

  cli::cli_h2("2) Scanning for rows with NA in critical columns")
  tbl <- remove_rows_with_na(tbl, columns, filepath = filepaths$removed_missing_critical)


  # 3) Remove duplicated SNPs --------------------------------------------------
  if(allow_duplications) {
    cli::cli_alert_warning("Skipping the check for duplicated variants, as `allow_duplications` is set to TRUE")

  } else {
    cli::cli_h2("3) Scanning for rows with duplications")
    tbl <- remove_duplicates(tbl, filepath = filepaths$removed_duplicates)

  }



  # MAP TO DBSNP ------------------------------------------------------------
  # -------------------------------------------------------------------------
  # -------------------------------------------------------------------------
  chr_pos_available <- all(c("CHR", "POS") %in% colnames(tbl))
  if(chr_pos_available) {

    # IF CHR and POS exists, we use that regardless of whether RSID is present
    # -------------------------------------------------------------------------
    cli::cli_h3("4) Validating columns")
    main_callback <- make_callback(filepaths$removed_validate_chr_pos)
    tbl <- validate_sumstat(tbl, filter_func = main_callback, convert_p = convert_p)

    cli::cli_h3("5) Adding RSID based on CHR:POS. Adding dbSNP based QC flags")
    if(build == "NA") {
      inferred_build <- infer_build(tbl, dbsnp_path = dbsnp_path)
    } else {
      inferred_build <- build
    }

    # -------------------------------------------------------------------------
    main <- repair_ids(tbl, build = inferred_build, dbsnp_path = dbsnp_path)

  } else {

    # no notion of build, as no CHR:POS
    inferred_build <- default_build

    # drop edge case of only CHR or only POS available
    tbl <- dplyr::select(tbl, -dplyr::any_of(c("CHR", "POS")))

    # update the RSID column
    tbl <- update_rsid(tbl, dbsnp_path = dbsnp_path, filepath = filepaths$updated_rsid)

    # check for CHR:POS in RSID column ----------------------------------------
    tbl <- validate_rsid(tbl, filepath = filepaths$removed_invalid_rsid)
    without_rsid <- tbl$without_rsid
    tbl <- tbl$main

    # validate columns --------------------------------------------------------

    cli::cli_h3("4a) Validating columns with a valid RSID")
    main_callback <- make_callback(filepaths$removed_validate_rsid)
    tbl <- validate_sumstat(tbl, filter_func = main_callback, convert_p = convert_p)

    cli::cli_h3("4b) Validating columns with CHR:POS in RSID column")
    without_rsid_callback <- make_callback(filepaths$removed_validate_rsid_without_rsid)
    without_rsid <- validate_sumstat(without_rsid, filter_func = main_callback, convert_p = convert_p)

    cli::cli_h3("5) Adding CHR and POS based on RSID. Adding dbSNP based QC flags")
    main <- repair_ids(tbl, dbsnp_path = dbsnp_path, repair = "pos")
    without_rsid <- repair_ids(without_rsid, dbsnp_path = dbsnp_path, repair = "rsid")
    main <- dplyr::bind_rows(main, without_rsid)


    # -------------------------------------------------------------------------
    # possible that rows coded as CHR:POS in the RSID column are actually
    # duplications of the other columns, but we cannot detect that until now
    if(allow_duplications) {
      cli::cli_alert_warning(
      "Skipping the check for duplicated variants.
      This means that variants with CHR:POS in RSID columns are allowed to have the same RSID as existing RSIDs in the dataset.
      Change this by setting `allow_duplications` to FALSE"
      )
    } else {
      cli::cli_h3("6) Scanning for rows with duplications")
      main <- remove_duplicates(
        main,
        columns = c("CHR" ,"POS_38", "EffectAllele","OtherAllele"),
        filepath = filepaths$removed_duplications_chr_pos_in_rsid_col
      )
    }


  }


  # DBSNP MAPPING DONE ------------------------------------------------------
  # -------------------------------------------------------------------------
  # -------------------------------------------------------------------------



  # apply filters -----------------------------------------------------------

  # no dbSNP mapping
  n_before <- nrow(main)
  before_filters <- main
  main <- dplyr::filter(main, !no_dbsnp_entry & !incompat_alleles) |>
    dplyr::select(-dplyr::all_of(c("no_dbsnp_entry", "incompat_alleles")))

  removed_no_dbsnp <- dplyr::anti_join(before_filters, main, by = "rowid")
  if(nrow(removed_no_dbsnp) > 0) {

    cli::cli_alert_warning("Removed {nrow(removed_no_dbsnp)} rows with no dbSNP entry or with incompat alleles")
    cli::cli_inform("{.file {filepaths$removed_no_dbsnp}}")
    arrow::write_parquet(removed_no_dbsnp, filepaths$removed_no_dbsnp)

  }

  main <-
    flag_duplicates(main, column = "rsid") |>
    dplyr::rename(multi_allelic = dup_rsid)


  # merge back indels -------------------------------------------------------

  if(!is.null(indels)) {
    if("POS" %in% colnames(indels)) {
      renamed_pos <- list("POS")
      names(renamed_pos) <- paste0("POS", "_", inferred_build)
      indels <- dplyr::rename(indels, !!!renamed_pos)
    }

    main <- dplyr::bind_rows(main, dplyr::mutate(indels, indel = TRUE))

  }


  # 6) repair missing statistics columns ------------------------------------

  if(repair_cols) {

    cli::cli_h2("6) Repairing missings statistics columns if possible")
    main <- repair_stats(main, impute_freq = impute_freq, impute_n = impute_n)

  }


  # all removed rows should be able to be tracked
  identify_removed_rows(dplyr::select(main,rowid), filepaths)


  # end of pipeline  ---------------------------------------------------------
  # -------------------------------------------------------------------------
  # -------------------------------------------------------------------------




  write_finished_tidyGWAS(df = main, output_format = output_format, filepaths = filepaths)
  fmt <- prettyunits::pretty_dt(Sys.time() - start_time)
  cli::cli_h1("Finished tidyGWAS")
  cli::cli_alert_info("A total of {rows_start - nrow(main)} rows were removed")
  cli::cli_alert_info("Total running time: {fmt}")

  metadata <- list(
    date = as.character(Sys.Date()),
    tidyGWAS_version = as.character(utils::packageVersion('tidyGWAS')),
    rows_start = rows_start,
    rows_end = nrow(main),
    basename = filename,
    raw_md5 = md5,
    dbsnp_path = dbsnp_path,
    output_format = output_format,
    build = build,
    outdir = output_dir,
    impute_freq = impute_freq,
    impute_n = impute_n,
    convert_p = convert_p,
    indel_strategy = indel_strategy,
    repair_cols = repair_cols,
    logfile = logfile,
    default_build = default_build,
    inferred_build = inferred_build,
    study_ControlN = ControlN,
    study_N = N,
    study_CaseN = CaseN
  )



  cli::cli_inform("Saving metadata from analysis to {.file {filepaths$metadata}}")
  metadata <- if(missing(column_names)) metadata else c(metadata, column_names)
  yaml::write_yaml(metadata, filepaths$metadata)

  # let function return the cleaned sumstats
  main

}




setup_pipeline_paths <- function(outdir, filename) {


  if(missing(filename)) filename <- "raw"
  stopifnot("The provided output folder already exists" = !dir.exists(outdir))
  pipeline_info <- paste(outdir,"pipeline_info", sep = "/")
  dir.create(outdir, recursive = TRUE)
  dir.create(pipeline_info, recursive = TRUE)
  dir.create(paste(outdir,"raw", sep = "/"))




  list(
    "base" = outdir,
    "logfile" = paste0(outdir, "/tidyGWAS_logfile.txt"),
    "cleaned" = paste(outdir, "tidyGWAS_hivestyle", sep = "/"),
    "raw_sumstats" = paste0(outdir, "/raw/", filename, ".parquet"),
    "metadata" = paste0(outdir, "/metadata.yaml"),
    "updated_rsid"= paste(pipeline_info, "updated_rsid.parquet", sep = "/"),

    "removed_missing_alleles" =                  paste(pipeline_info, "removed_missing_alleles.parquet", sep = "/"),
    "removed_indels" =                           paste(pipeline_info, "removed_indels.parquet", sep = "/"),
    "removed_validate_indels" =                  paste(pipeline_info, "removed_validate_indels.parquet", sep = "/"),

    "removed_missing_critical" =                 paste(pipeline_info, "removed_missing_critical.parquet", sep = "/"),
    "removed_duplicates" =                       paste(pipeline_info, "removed_duplicates.parquet", sep = "/"),

    "removed_invalid_rsid" =                     paste(pipeline_info, "removed_invalid_rsid.parquet", sep = "/"),
    "removed_validate_rsid" =                    paste(pipeline_info, "removed_validate_rsid_path.parquet", sep = "/"),
    "removed_validate_rsid_without_rsid" =       paste(pipeline_info, "removed_without_rsid.parquet", sep ="/"),
    "removed_duplications_chr_pos_in_rsid_col" = paste(pipeline_info, "removed_duplications_chr_pos_in_rsid_col.parquet", sep = "/"),

    "removed_validate_chr_pos" =                 paste(pipeline_info, "removed_validate_chr_pos_path.parquet", sep = "/"),

    "removed_no_dbsnp" =                         paste(pipeline_info, "removed_nodbsnp.parquet", sep = "/")
  )


}
