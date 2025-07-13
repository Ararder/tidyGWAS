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
#' @param dbsnp_path filepath to the dbSNP155 directory
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
#' @param impute_freq one of c("None", "EUR", "AMR", "AFR", "SAS", "EAS"). If None, no imputation is done.
#'  Otherwise precomputed alleles frequence from 1000KG, selected ancestry is used
#' @param impute_freq_file filepath to a .parquet file with custom allele frequencies.
#'  The file needs to be a tabular dataframe with columns RSID, EffectAllele, OtherAllele, EAF.
#'  EAF should correspond to the frequency of the EffectAllele.
#' @param min_EAF Apply a filter on allele frequency prior to applying the algorithm. Useful to speed up cleaning of very large files
#' @param flag_discrep_freq Should variants with allele frequency discrepancies be flagged?
#' @param impute_n Should N be imputed if it's missing?
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
#' tidyGWAS(
#'   tbl = "path/to/GWAS_trait_X_.tsv.gz", logfile = TRUE,
#'   output_dir = "/store/GWAS/tidyGWAS/trait_X"
#'   )
#' }
tidyGWAS <- function(
    tbl,
    dbsnp_path,
    ...,
    column_names = NULL,
    output_format = c("hivestyle","parquet", "csv"),
    output_dir = tempfile(),
    CaseN = NULL,
    ControlN = NULL,
    N = NULL,
    impute_freq = c("None", "EUR", "AMR", "AFR", "SAS", "EAS"),
    impute_freq_file = NULL,
    impute_n = FALSE,
    min_EAF = NULL,
    flag_discrep_freq = c("None", "EUR", "AMR", "AFR", "SAS", "EAS"),
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
  (rlang::is_scalar_character(tbl) | "data.frame" %in% class(tbl)) || cli::cli_abort(
    "tbl must be a character vector with length 1, or a data.frame"
  )


  rlang::check_required(dbsnp_path)
  rlang::is_scalar_character(dbsnp_path) || cli::cli_abort("dbsnp_path must be a single character string")
  file.exists(dbsnp_path) || cli::cli_abort("The dbsnp_path provided does not exist")


  if(!is.null(CaseN)) stopifnot(rlang::is_scalar_integerish(CaseN))
  if(!is.null(ControlN))  stopifnot(rlang::is_scalar_integerish(ControlN))
  if(!is.null(N)) stopifnot(rlang::is_scalar_integerish(N))
  impute_freq <- rlang::arg_match(impute_freq)
  flag_discrep_freq <- rlang::arg_match(flag_discrep_freq)
  if(!is.null(impute_freq_file) & impute_freq != "None")  {
    cli::cli_abort(
      "If `impute_freq_file` is provided, `impute_freq` must be `None`.
      Impute_freq is used to add allele frequency from the provided 1000kg reference files.
      impute_freq_file is used to pass a custom file - Using both are imcompatible."
    )
  }
  if(!is.null(impute_freq_file)) stopifnot(rlang::is_scalar_character(impute_freq) && file.exists(impute_freq))


  if(!is.null(min_EAF)) {
    rlang::is_scalar_double(min_EAF) || rlang::abort("min_EAF must be a single value of type double")
    if(min_EAF <= 0 | min_EAF >= 0.5) rlang::abort("min_EAF must be a double in the range 0 <= min_EAF <= 0.5")

  }

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


  # -------------------------------------------------------------------------


  # welcome message ----------------------------------------------------------
  start_time <- Sys.time()

  # check input data.frame --------------------------------------------------
  cli::cli_inform("Parsing input summary statistics...")
  tbl <- parse_input(tbl, ...)
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
  # welcome message ----------------------------------------------------------
  start_time <- Sys.time()
  cli::cli_h1("Running {.pkg tidyGWAS {packageVersion('tidyGWAS')}}")
  cli::cli_inform("Starting at {start_time}")
  cli::cli_inform("with {rows_start} rows in input data.frame")
  cli::cli_alert_info("Saving output in folder: {.file {filepaths$base}}")

  # The sumstats are saved without any edits, to not loose information
  arrow::write_parquet(tbl,  filepaths$raw_sumstats)


  # start of pipeline ----------------------------------------------------------
  if(!is.null(CaseN)) tbl <- dplyr::mutate(tbl, CaseN =  {{ CaseN }})
  if(!is.null(ControlN)) tbl <- dplyr::mutate(tbl, ControlN =  {{ ControlN }})
  if(!is.null(N)) tbl <- dplyr::mutate(tbl, N = {{ N }})

  # 0) formatting --------------------------------------------------------------
  if(is.null(column_names)) {
    tbl <- guess_names(tbl)

  } else {
    check_columns(column_names, tbl)
    tbl <- dplyr::rename(tbl, !!!column_names)
  }


  tbl <- select_correct_columns(tbl)

  if(!is.null(min_EAF)) {
    "EAF" %in% colnames(tbl) || rlang::abort("EAF column is missing from the data.frame. Cannot prefilter on allele frequency without EAF")
    cli::cli_inform("Filtering rows with EAF < {min_EAF}")
    n_eaf <- nrow(tbl)
    tbl <- dplyr::filter(tbl, EAF >= min_EAF & EAF <= (1-min_EAF))
    n_eaf_after <- nrow(tbl)
    cli::cli_alert_info("Removed {n_eaf - n_eaf_after} rows with EAF < {min_EAF}")
  }



  # 1) detect indels ----------------------------------------------------------
  cli::cli_h2("1) Scanning for indels")

  tbl <- remove_rows_with_na(tbl, c("EffectAllele", "OtherAllele"), filepaths$removed_missing_alleles)
  tbl <- detect_indels(tbl, indel_strategy = indel_strategy, filepaths = filepaths, convert_p = convert_p)

  indels <- tbl$indels
  if(nrow(indels) == 0) indels <- NULL
  tbl <- tbl[["main"]]


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
    main <- apply_dbsnp_filter(main, filepaths)

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
    cli::cli_inform("Number of rows with valid RSID: {nrow(tbl)}")
    main_callback <- make_callback(filepaths$removed_validate_rsid)
    tbl <- validate_sumstat(tbl, filter_func = main_callback, convert_p = convert_p)
    main <- repair_ids(tbl, dbsnp_path = dbsnp_path, repair = "pos")

    if(!is.null(without_rsid)) {
      cli::cli_h3("4b) Validating columns with CHR:POS in RSID column")
      cli::cli_inform("Number of rows with CHR:POS in RSID column: {nrow(without_rsid)}")

      without_rsid_callback <- make_callback(filepaths$removed_validate_rsid_without_rsid)
      without_rsid <- validate_sumstat(without_rsid, filter_func = main_callback, convert_p = convert_p)
      without_rsid <- repair_ids(without_rsid, dbsnp_path = dbsnp_path, repair = "rsid")

      main <- dplyr::bind_rows(main, without_rsid)
    }

    main <- apply_dbsnp_filter(main, filepaths)

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

    main <- dplyr::bind_rows(
      dplyr::mutate(main, indel = FALSE),
      dplyr::mutate(indels, indel = TRUE)
    )

  }


  # 6) repair missing statistics columns ------------------------------------

  if(repair_cols) {

    cli::cli_h2("6) Repairing missings statistics columns if possible")
    main <- repair_stats(main,dbsnp_path = dbsnp_path, impute_freq = impute_freq, impute_freq_file = impute_freq_file, impute_n = impute_n)

  }


  # all removed rows should be able to be tracked
  identify_removed_rows(dplyr::select(main,rowid), filepaths)


  if(flag_discrep_freq != "None") {
    main <- add_freq_diff_flag(main, flag_discrep_freq,dbsnp_path)

  }


  # end of pipeline  ---------------------------------------------------------
  # -------------------------------------------------------------------------
  # -------------------------------------------------------------------------




  write_finished_tidyGWAS(df = main, output_format = output_format, filepaths = filepaths)
  fmt <- prettyunits::pretty_dt(Sys.time() - start_time)
  cli::cli_h1("Finished tidyGWAS")
  cli::cli_alert_info("A total of {rows_start - nrow(main)} rows were removed")
  cli::cli_alert_info("Total running time: {fmt}")
  # args <- modifyList(formals(), as.list(match.call()))
  # implement later
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
