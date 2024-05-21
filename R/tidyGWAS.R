utils::globalVariables(c(
  "chr_pos", "n", "n_chr_pos", "dup_rsid",
  ".data", "invalid_rsid",  "new_rsid",
  "has_rsid", "head", "CHR_37", "POS_37",
   "reason", "filter_callback", "no_dbsnp_entry", "logfile",
  "alt_allele"))

snp_cols <- c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele", "rowid")
info_cols <- c("INFO", "N", "CaseN", "ControlN", "EAF")
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
#' the correct type. Missing CHR/POS or RSID is detected and imputed using
#' [repair_rsid()] or [repair_chr_pos()].
#'
#' If statistis such as `P`, `B` are missing, `tidyGWAS()` will attempt to impute
#' them if possible using [repair_stats()]
#'
#' Standard column names are assumed, BEFORE inputting into the function. This is
#' a deliberate decision, as automatic parsing of some important column names
#' can be ambigious For example, in some sumstats, A1 referes to effect allele,
#' while other formats use A1 as non-effect allele. [tidyGWAS_columns()] can be
#' used to standardise column names, and see the standard format.
#'
#' @param tbl a `data.frame` or `character()` vector
#' @param dbsnp_path filepath to the dbSNP155 directory (untarred dbSNP155.tar)
#' @param ... pass additional arguments to [arrow::read_delim_arrow()], if tbl is a filepath.
#' @param output_dir filepath to a folder where data should be stored. The folder should not yet exist.
#' @param output_format How should the finished cleaned file be saved?
#'  * "'csv' corresponds to [arrow::write_csv_arrow()]
#'  * 'parquet' corresponds to [arrow::write_parquet()]
#'  * 'hivestyle' corresponds to [arrow::write_dataset()] split by CHR
#'
#' @param logfile Should messages be redirected to a logfile?
#' @param column_names a named list of column names, to be used to rename columns.
#' @param CaseN manually input number of cases
#' @param ControlN manually input number of controls
#' @param N manually input sample size
#' @param build If you are sure of what genome build ('37' or '38'), can be used to skip [infer_build()] and speed up computation
#' @param convert_p What value should be used for when P-value has been rounded to 0?
#' @param indel_strategy Should indels be kept or removed?
#' @param repair_cols Should any missing statistical columns be repaired if possible?
#' @param add_missing_build Should the build which the sumstats are NOT on also be added in?
#' @param verbose Explain filters in detail?
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
    output_dir = paste0(tempdir(), "/",stringr::str_replace_all(date(), pattern = c(" "="_", ":"="_"))),
    CaseN = NULL,
    ControlN = NULL,
    N = NULL,
    build = c("NA","37", "38"),
    indel_strategy = c("keep", "remove"),
    convert_p = 2.225074e-308,
    repair_cols = TRUE,
    logfile = FALSE,
    verbose = FALSE,
    add_missing_build = TRUE
    ) {



  # parse arguments ---------------------------------------------------------
  rlang::check_required(tbl)
  rlang::check_required(dbsnp_path)
  if(!is.null(CaseN)) stopifnot(rlang::is_scalar_integerish(CaseN))
  if(!is.null(ControlN))  stopifnot(rlang::is_scalar_integerish(ControlN))
  if(!is.null(N)) stopifnot(rlang::is_scalar_integerish(N))
  stopifnot(rlang::is_scalar_double(convert_p))
  stopifnot(rlang::is_bool(repair_cols))
  stopifnot(rlang::is_bool(logfile))
  stopifnot(rlang::is_bool(verbose))
  stopifnot(rlang::is_bool(add_missing_build))
  stopifnot("logfile can only be TRUE or FALSE"= rlang::is_bool(logfile))
  stopifnot("The `output_dir` specified already exists" = !dir.exists(output_dir))
  stopifnot("the filepath for dbSNP does not exist" = file.exists(dbsnp_path))
  start_time <- Sys.time()
  output_format <-  rlang::arg_match(output_format)
  build <-          rlang::arg_match(build)
  indel_strategy <- rlang::arg_match(indel_strategy)


  # welcome message ----------------------------------------------------------
  cli::cli_h1("Running {.pkg tidyGWAS {packageVersion('tidyGWAS')}}")
  cli::cli_inform("Starting at {start_time}")






  # check input data.frame --------------------------------------------------

  tbl <- parse_tbl(tbl, ...)
  filename <- tbl$filename
  md5 <- tbl$md5
  # has to be last, overwriting tbl here
  tbl <- tbl$tbl

  rows_start <- nrow(tbl)
  filepaths <- setup_pipeline_paths(outdir = output_dir, filename = filename)

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


  # 1) Drop na -----------------------------------------------------------------

  cli::cli_h2("1) Scanning for rows with NA")
  tbl <- remove_rows_with_na(tbl, filepaths = filepaths)


  # 2) Remove duplicated SNPs --------------------------------------------------

  cli::cli_h2("2) Scanning for rows with duplications")
  tbl <- remove_duplicates(tbl, filepaths = filepaths)

  # 3) detect indels ----------------------------------------------------------

  cli::cli_h2("3) Scanning for indels")
  tbl <- detect_indels(tbl, indel_strategy = indel_strategy, filepaths = filepaths, verbose = verbose, convert_p = convert_p)

  indels <- tbl$indels
  tbl <- tbl$main


  # 4) update/repair RSID ---------------------------------------------------


  # if only RSID exists, several steps needs to be taken, see rsid_only()
  if("RSID" %in% colnames(tbl) & !all(c("CHR", "POS") %in% colnames(tbl))) {

    main <- rsid_only(
      tbl,
      dbsnp_path = dbsnp_path,
      filepaths = filepaths,
      verbose = verbose,
      convert_p = convert_p,
      add_missing_build = add_missing_build
    )
    # if only RSID is present, there is not a notion of 'inferred_build',
    # since RSID is supposed to be build agnostic
    inferred_build <- NULL


  # if CHR:POS exists, we use that, regardless of whether RSID exists
  } else {

    cli::cli_h3("4a) Validating columns")
    main_callback <- make_callback(filepaths$removed_validate_chr_pos)
    tbl <- validate_sumstat(tbl, filter_func = main_callback, verbose = verbose, convert_p = convert_p)


    # -------------------------------------------------------------------------

    cli::cli_h3("5) Adding RSID based on CHR:POS. Adding dbSNP based QC flags")

    if(build == "NA") {
      inferred_build <- infer_build(tbl, dbsnp_path = dbsnp_path)
    } else {
      inferred_build <- build
    }
    main <- repair_rsid(tbl, build = inferred_build, dbsnp_path = dbsnp_path, add_missing_build = add_missing_build)

  }


  # apply filters -----------------------------------------------------------

  main <- apply_filters(tbl = main, filepaths = filepaths) |>
    flag_duplicates(column = "rsid") |>
    dplyr::rename(multi_allelic = dup_rsid)


  # merge back indels -------------------------------------------------------

  if(!is.null(indels)) {

    main <- dplyr::bind_rows(main, dplyr::mutate(indels, indel = TRUE))

  }


  # 6) repair missing statistics columns ------------------------------------

  if(repair_cols) {

    cli::cli_h2("6) Repairing missings statistics columns if possible")
    main <- repair_stats(main)

  }

  # add a CHR:POS:REF:ALT column
  main <- create_id(main)

  # all removed rows should be able to be tracked
  identify_removed_rows(dplyr::select(main,rowid), filepaths)


  # end of pipeline  ---------------------------------------------------------



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
    convert_p = convert_p,
    indel_strategy = indel_strategy,
    repair_cols = repair_cols,
    logfile = logfile,
    verbose = verbose,
    add_missing_build = add_missing_build,
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





# -------------------------------------------------------------------------
# -------------------------------------------------------------------------




#' Parse the input data.frame or filepath to tidyGWAS
#' @inheritParams tidyGWAS
#' @param ... optional arguments passed to [arrow::read_delim_arrow()]
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' # read in file from disk, using tab delimiter and skipping first 100 rows
#' df <- parse_tbl("path/to_sumstats", delim = "\t", skip = 100)
#' df <- parse_tbl(tibble(RSID = "rs1001", B = 0.01, P = 0.005, SE = 0.001, N = 100))
#' }
parse_tbl <- function(tbl, ...) {

  rlang::check_required(tbl)


  if("character" %in% class(tbl)) {
    stopifnot("A filepath has to be a character vector of length 1" = "character" %in% class(tbl) & length(tbl) == 1)
    stopifnot("File does not exist"  = file.exists(tbl))
    filename <- basename(tbl)
    md5 <- tools::md5sum(tbl)
    tbl <- arrow::read_delim_arrow(tbl, ...)

    ncols <- length(colnames(tbl))
    if(ncols <= 3) {
      stop(glue::glue("Only {ncols} columns in the file. At least 3 columns are required.
                      Most likely you forgot to use the `delim` argument to specify the delimiter."))
    }


  } else if("data.frame" %in% class(tbl)) {
    md5 <- NULL
    tbl <- dplyr::tibble(tbl)
    filename <- "raw"

  } else {

    stop("tbl is not a character vector or a data.frame")
  }

  if(!"rowid" %in% colnames(tbl)) tbl$rowid <- 1:nrow(tbl)

  list("tbl" = tbl, "filename" = filename, "md5" = md5)
}





# -------------------------------------------------------------------------
# -------------------------------------------------------------------------




identify_removed_rows <- function(finished, filepaths) {

  start <- arrow::read_parquet(filepaths$raw_sumstats, col_select = "rowid")

  removed_rows <- dplyr::anti_join(start, finished, by = "rowid")

  # remove the file with all rowids
  files_in_dir <- list.files(paste0(filepaths$base, "/pipeline_info/"), pattern = "*removed_*", full.names = TRUE)


  removed_rows_with_flags <-
    files_in_dir |>
    purrr::set_names(base::basename) |>
    purrr::map(arrow::read_parquet) |>
    purrr::map(\(x) dplyr::select(x, rowid)) |>
    purrr::list_rbind(names_to = "reason") |>
    dplyr::select(rowid, reason) |>
    dplyr::mutate(reason = stringr::str_remove(reason, "removed_")) |>
    dplyr::mutate(reason = stringr::str_remove(reason, ".parquet"))

  if(sum(!removed_rows$rowid %in% removed_rows_with_flags$rowid)) {
    cli::cli_alert_danger(
      "WARNING: Could not track why some rows were removed. This means that some rows were removed
      that has not been accounted for. possibly unintended rows were removed.")
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



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------



#' Create the folder structure for tidyGWAS
#'
#' @param outdir output directory
#' @param filename filename for the raw sumstats
#' @return a list of filepaths
#' @export
#'
#' @examples
#' setup_pipeline_paths(tempfile())
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

    "removed_rows"=                              paste(pipeline_info, "removed_", sep = "/"),
    "removed_duplicates" =                       paste(pipeline_info, "removed_duplicates.parquet", sep = "/"),
    "failed_rsid_parse" =                        paste(pipeline_info, "removed_failed_rsid_parse.parquet", sep = "/"),
    "removed_invalid_chr_pos_in_rsid"  =         paste(pipeline_info, "removed_invalid_chr_pos_rsid.parquet", sep = "/"),
    "removed_indels" =                           paste(pipeline_info, "removed_indels.parquet", sep = "/"),
    "removed_validate_rsid" =                    paste(pipeline_info, "removed_validate_rsid_path.parquet", sep = "/"),
    "removed_validate_rsid_without_rsid" =       paste(pipeline_info, "removed_without_rsid.parquet", sep ="/"),
    "removed_validate_chr_pos" =                 paste(pipeline_info, "removed_validate_chr_pos_path.parquet", sep = "/"),
    "removed_validate_indels" =                  paste(pipeline_info, "removed_validate_indels.parquet", sep = "/"),
    "removed_duplications_chr_pos_in_rsid_col" = paste(pipeline_info, "removed_duplications_chr_pos_in_rsid_col.parquet", sep = "/"),
    "removed_no_dbsnp" =                         paste(pipeline_info, "removed_nodbsnp.parquet", sep = "/"),
    "removed_chr_mismatch" =                     paste(pipeline_info, "removed_chr_mismatch.parquet", sep = "/"),
    "removed_missing_on_either_build" =          paste(pipeline_info, "removed_missing_on_either_build.parquet", sep = "/")
  )





}


write_finished_tidyGWAS <- function(df, output_format, filepaths) {

  df <- standardize_column_order(df)

  if(output_format == "hivestyle") {

    arrow::write_dataset(dplyr::group_by(df, CHR), filepaths$cleaned)

  } else if(output_format == "parquet") {

    arrow::write_parquet(df, paste(filepaths$base, "tidyGWAS_cleaned.parquet", sep = "/"))

  } else if(output_format == "csv") {

    arrow::write_csv_arrow(df, paste(filepaths$base, "tidyGWAS_cleaned.csv", sep = "/"))

  }


}



# -------------------------------------------------------------------------


#' Create a [dplyr::tibble()] with tidyGWAS column names
#'
#' [tidyGWAS()] requires the column names to be in a specific format.
#' This function facilitates that format, and also provides a list of all the columns
#' tidyGWAS considers valid. Any columns with a non tidyGWAS column name will be dropped.
#'
#' The function is a simple wrapper around [dplyr::select()], with each possible
#' tidyGWAS column as an argument-
#'
#' @param tbl a [dplyr::tibble()] or something coercible to one
#' @param CHR chromosome
#' @param POS position
#' @param RSID rsID from dbSNP
#' @param EffectAllele allele corresponding to effect
#' @param OtherAllele  non-effect allele
#' @param B Beta, effect,
#' @param SE standard error
#' @param P p value
#' @param EAF effect-allele frequency
#' @param N total sample size (case + control)
#' @param CaseN number of cases (for case-control phenotypes)
#' @param ControlN number of controls (for case-control phenotypes)
#' @param INFO INFO score, imputation accuracy
#' @param Z Z score
#' @param OR odds-ratio
#'
#' @return a tibble with changes column names
#' @export
#'
#' @examples
#' wrong_format <- dplyr::tibble(CHROM = 1, bp = 1000, A1 = "C", A2 = "A", Effect = 0.05)
#' formatted <- tidyGWAS_columns(
#' wrong_format, CHR = "CHROM", POS = "bp",
#' EffectAllele = "A1", OtherAllele = "A2", B = "Effect"
#' )
#' # columns that are wrongly named and NOT passed to tidyGWAS_columns() are dropped
#' wrong_format <- dplyr::tibble(CHROM = 1, bp = 1000, A1 = "C", A2 = "A", Effect = 0.05)
#' tidyGWAS_columns(
#' wrong_format, CHR = "CHROM", POS = "bp",
#' EffectAllele = "A1"
#' )
#'
#'
tidyGWAS_columns <- function(
    tbl,
    CHR = "CHR",
    POS = "POS",
    RSID = "RSID",
    EffectAllele = "EffectAllele",
    OtherAllele = "OtherAllele",
    B = "B",
    SE = "SE",
    P = "P",
    EAF = "EAF",
    N = "N",
    CaseN = "CaseN",
    ControlN = "ControlN",
    INFO = "INFO",
    Z = "Z",
    OR = "OR"
) {



  tbl |>
    dplyr::select(dplyr::any_of(c(
      "CHR" = CHR,
      "POS" = POS,
      "RSID" = RSID,
      "EffectAllele" = EffectAllele,
      "OtherAllele" = OtherAllele,
      "B" = B,
      "SE" = SE,
      "P" = P,
      "EAF" = EAF,
      "N" = N,
      "CaseN" = CaseN,
      "ControlN" = ControlN,
      "INFO" = INFO,
      "Z" = Z,
      "OR" = OR
    )))

}


standardize_column_order <- function(tbl) {
  dplyr::select(tbl, dplyr::any_of(
    c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele", "EAF",
      "Z", "B", "SE", "P", "N", "CaseN", "ControlN", "INFO"
      ,"CHR_37", "POS_37", "rowid", "multi_allelic", "indel", "REF" = "ref_allele")), dplyr::everything()
    )

}

check_correct_files <- function(dir) {


  stopifnot(
    "The reference files have been damaged or are missing. Please re-download the reference files.,
        Expected to find folders GRCh37, GRCh38 and refsnp-merged in reference dir." =
      all(c("GRCh37", "GRCh38", "refsnp-merged") %in% list.files(dir))
  )

  stopifnot(
    "Expected to find a total of 51 parquet files across GRCh37 and GRCh38 and merged RSIDs. This
    indicates that the reference files have been damaged or corrupted. Please re-download the reference files" =
      length(list.files(dir, recursive = TRUE)) == 51
  )


  cli::cli_alert_success("Reference files are present and correct. Saved in {.path {dir}}")

}



rsid_only <- function(tbl, dbsnp_path, filepaths, verbose, convert_p, add_missing_build) {

  # update the RSID column
  tbl <- update_rsid(tbl, dbsnp_path = dbsnp_path, filepaths = filepaths)


  # check for CHR:POS in RSID column ----------------------------------------
  tbl <- validate_rsid(tbl, outpath = filepaths$removed_invalid_chr_pos_in_rsid)
  without_rsid <- tbl$without_rsid
  tbl <- tbl$main


  # validate columns --------------------------------------------------------
  cli::cli_h3("4a) Validating columns with a correct RSID")

  main_callback <- make_callback(filepaths$removed_validate_rsid)
  tbl <- validate_sumstat(tbl, filter_func = main_callback, verbose = verbose, convert_p = convert_p)


  # We have the validate the data.frame we split of, that had CHR:POS in RSID
  cli::cli_h3("4b) Validating columns with CHR:POS in RSID column")
  without_rsid_callback <- make_callback(filepaths$removed_validate_rsid_without_rsid)
  without_rsid <- validate_sumstat(without_rsid, filter_func = main_callback, verbose = verbose, convert_p = convert_p)


  cli::cli_h3("5) Adding CHR and POS based on RSID. Adding dbSNP based QC flags")

  main <- repair_chr_pos(tbl, dbsnp_path = dbsnp_path, add_missing_build = add_missing_build)
  without_rsid <- repair_rsid(without_rsid, dbsnp_path = dbsnp_path, add_missing_build = add_missing_build)

  # merge the data.frames together
  main <- dplyr::bind_rows(main, without_rsid)


  # check for dups ----------------------------------------------------------
  # It's possible that rows coded as CHR:POS in the RSID column are actually
  # duplications of the other columns, but we cannot detect that until we
  # have CHR:POS:RSID for all variants. So we have to do a check here.


  before_unique_check <- dplyr::select(main, rowid)

  # duplication check fails when CHR:POS is missing, so have split by builds
  b38_missing <- dplyr::filter(main, is.na(CHR)) |>
    dplyr::distinct(CHR_37,POS_37,EffectAllele,OtherAllele, .keep_all = TRUE)

  main <- dplyr::filter(main, !is.na(CHR)) |>
    dplyr::distinct(CHR,POS,EffectAllele,OtherAllele, .keep_all = TRUE) |>
    dplyr::bind_rows(b38_missing)


  removed <- dplyr::anti_join(before_unique_check, main, by = "rowid")

  if(nrow(removed) > 0) {
    cli::cli_alert_info(
    "Found {nrow(removed)} rows which are duplicates. These duplicates comes from the subset of
    rows that had CHR:POS in the RSID column. They could not be detected as duplicates before
    their CHR and POS was inferred.")
    arrow::write_parquet(removed, filepaths$removed_duplications_chr_pos_in_rsid_col)
  }

  # -------------------------------------------------------------------------

  main

}

update_column_names <- function(tbl, column_map, CaseN = NULL, ControlN =NULL, N=NULL) {
  rlang::check_required(tbl)


  if(!missing(column_map)) {


    stopifnot(
        "All entries in column_map must be present in the input tbl" =
          all(column_map %in% colnames(tbl))
        )

    tbl <- dplyr::rename(tbl, !!!column_map)
  }

  if(!is.null(CaseN)) tbl <- dplyr::mutate(tbl, CaseN =  {{ CaseN }})
  if(!is.null(ControlN)) tbl <- dplyr::mutate(tbl, ControlN =  {{ ControlN }})
  if(!is.null(N)) tbl <- dplyr::mutate(tbl, N = {{ N }})





  tbl
}



check_zero_rows <- function(tbl){
  if(!is.null(tbl)) {
    if(nrow(tbl) == 0) {
      tbl <- NULL
    }
  }

  return(tbl)
}

create_id <- function(tbl) {
  tbl |>
    dplyr::mutate(
      alt_allele = dplyr::if_else(ref_allele == EffectAllele, OtherAllele, EffectAllele),
      ID = stringr::str_c(CHR, POS, ref_allele, alt_allele, sep = ":")
    ) |>
    dplyr::select(-"alt_allele")

}

