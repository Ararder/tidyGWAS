utils::globalVariables(c(
  "chr_pos", "n", "n_chr_pos", "dup_rsid",
  ".data", "invalid_rsid",  "new_rsid",
  "has_rsid", "head", "CHR_37", "POS_37",
   "reason", "filter_callback", "no_dbsnp_entry", "logfile"))

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
#' [repair_rsid()] or [repair_chr_pos()]. If both RSID and CHR:POS is present,
#' [verify_chr_pos_rsid()] is executed to check that dbSNP CHR:POS:RSID agrees
#' with CHR:POS:RSID in `tbl`.
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
#' @param ... pass additional arguments to [arrow::read_delim_arrow()], if tbl is a filepath.
#' @param output_format How should the finished cleaned file be saved?
#'  * "csv" corresponds to [arrow::write_csv_arrow()]
#'  * 'hivestyle' corresponds to [arrow::write_dataset()] split by CHR
#'  * 'parquet' corresponds to [arrow::write_parquet()]
#'
#' @param build Can be used to skip [infer_build()]
#' @param outdir Where should results be saved after a succesful execution?
#' @param convert_p What value should be used for P = 0?
#' @param name name of the output directory
#' @param dbsnp_path filepath to the dbSNP155 directory (untarred dbSNP155.tar)
#' @param study_n Sometimes N is missing from GWAS summary statistics. It is then
#' often much more useful to set a study-wide N for all rows, instead of leaving
#' the N column missing. study_n can be used to set the N column.
#' @param keep_indels Should indels be kept?
#' @param repair_cols Should any missing columns be repaired?
#' @param add_missing_build Should the build which the sumstats are NOT on also be added in?
#' @param logfile Should messages be redirected to a logfile?
#' @param log_on_err Optional. Can pass a filepath to copy the logfile to when the function exists.
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
    ...,
    dbsnp_path,
    output_format = c("csv","hivestyle", "parquet"),
    build = c("NA","37", "38"),
    outdir = tempdir(),
    study_n,
    convert_p = 2.225074e-308,
    name = stringr::str_replace_all(date(), pattern = c(" "="_", ":"="_")),
    keep_indels = TRUE,
    repair_cols = TRUE,
    logfile = FALSE,
    log_on_err="tidyGWAS_logfile.txt",
    verbose = FALSE,
    add_missing_build = TRUE
    ) {


  # parse arguments ---------------------------------------------------------
  start_time <- Sys.time()
  tbl <- parse_tbl(tbl, ...)
  rows_start <- nrow(tbl)
  output_format <- rlang::arg_match(output_format)
  build = rlang::arg_match(build)
  filepaths <- setup_pipeline_paths(name = name)

  # setup logging -----------------------------------------------------------
  if(isTRUE(logfile)) {
    cli::cli_alert_info("Output is redirected to logfile: {.file {filepaths$logfile}}")
    withr::local_message_sink(filepaths$logfile)
    withr::local_output_sink(filepaths$logfile)
    if(!missing(log_on_err)) on.exit(file.copy(filepaths$logfile, log_on_err), add=TRUE)
  }

  # write out the raw sumstats to always be able to find what changes was made to input file
  arrow::write_parquet(tbl, paste0(filepaths$base , "/raw_sumstats.parquet"), compression = "gzip")

  # welcome message ----------------------------------------------------------
  cli::cli_h1("Running {.pkg tidyGWAS {packageVersion('tidyGWAS')}}")
  cli::cli_inform("Starting at {start_time}, with {rows_start} rows in input data.frame")
  cli::cli_alert_info("Saving all files during execution to {.file {filepaths$base}},")
  cli::cli_alert_info("After execution, files will be copied to {.file {paste(outdir,name, sep = '/')}}")

  # remove columns without tidyGWAS column names
  tbl <- select_correct_columns(tbl, study_n)

  # drop any rows with NAs
  cli::cli_h2("1) Scanning for rows with NA")
  tbl <- remove_rows_with_na(tbl, filepaths)

  # update RSID
  cli::cli_h2("2) Updating merged RSIDs")
  if("RSID" %in% colnames(tbl) & !missing(dbsnp_path)) tbl <- update_rsid(tbl, dbsnp_path = dbsnp_path, filepaths = filepaths)

  # remove duplicated rows, using CHR:POS:REF:ALT or RSID:REF:ALT to construct ID
  cli::cli_h2("3) Scanning for rows with duplications")
  tbl <- remove_duplicates(tbl, filepaths = filepaths)


  # -------------------------------------------------------------------------
  # Need to split out the sumstats into separate data.frames
  # 1) main, 2) indels, 3) inalid_rsid

  data_list <- vector("list", length = 3)
  names(data_list) <- c("main", "indels", "without_rsid")

  cli::cli_h2("4) Scanning for indels")
  tbl <- detect_indels(tbl, keep_indels, filepaths)
  data_list$main <- tbl$main
  if(!is.null(tbl$indels)) data_list$indels <- tbl$indels
  tbl <- NULL

  # this function splits apart the data.frame into rows with and withot valid RSID
  # as rows without RSID will need to go into another dbSNP cleaning function
  if("RSID" %in% colnames(data_list$main)) {

    cli::cli_h2("5) Scanning for invalid RSID ")
    data_list$main <- validate_rsid(data_list$main, outpath = paste(filepaths$removed_rows, "invalid_chr_pos_rsid.parquet"))
    if(!is.null(data_list$main$without_rsid)) data_list$without_rsid <- data_list$main$without_rsid
    data_list$main <- data_list$main$main
    if(nrow(data_list$main) == 0) data_list[1] <- list(NULL)

  }

  # here the cleaning is split apart into three separate data.tables: main, indels, without_rsid
  cli::cli_h2("6) Column validation is done separately for main rows, rows without RSID and indels")
  filter_funcs <-  purrr::map(paste0(filepaths$removed_rows, c("main", "indels", "without_rsid"), "_col_validation"), make_callback)
  cols_to_not_validate <- list("", c("EffectAllele","OtherAllele"), "")
  id <- list("main rows", "indel rows", "rows without RSID")

  # validate the columns that it is possible to validate for each separate data.frame
  data_list <- list(tbl = data_list, remove_cols = cols_to_not_validate, filter_func = filter_funcs, verbose = list(verbose), convert_p = list(convert_p), id = id) |>
    purrr::pmap(validate_sumstat)

  # perform checks that need dbSNP if dbsnp_paths is passed
  if(!missing(dbsnp_path)) {

    cli::cli_h2("7 Using dbSNP to repair and validate CHR:POS:RSID")
    cli::cli_li("Variants with CHR and POS or RSID not present in dbSNP are removed")
    cli::cli_li("Checking that EffectAllele and OtherAllele is compatible with REF and ALT in dbSNP")

    data_list <- validate_with_dbsnp(data_list = data_list, build = build, dbsnp_path = dbsnp_path, filepaths = filepaths, add_missing_build = add_missing_build)

    main <- data_list

  } else {
    main <- dplyr::bind_rows(data_list$main, data_list$indels)
  }

  data_list <- NULL

  # repair stats if passed
  if(repair_cols) {

    cli::cli_h2("8) Repairing missings columns if possible")
    main <- repair_stats(main)

  }

  # flag indels and multi-allelics
  main <-
    flag_duplicates(main, "rsid") |>
    flag_indels() |>
    dplyr::rename(multi_allelic = dup_rsid)

  # make sure all removed rows can be tracked
  identify_removed_rows(dplyr::select(main,rowid), filepaths)


  # end of pipeline  --------------------------------------------------------
  # write out into fileformat dependnding on output_format
  write_finished_tidyGWAS(df = main, output_format = output_format, outdir = outdir, filepaths = filepaths)

  # exit header
  fmt <- prettyunits::pretty_dt(Sys.time() - start_time)
  cli::cli_h1("Finished tidyGWAS")
  cli::cli_alert_info("A total of {rows_start - nrow(main)} rows were removed")
  cli::cli_alert_info("Total running time: {fmt}")

  # -------------------------------------------------------------------------
  if(outdir != tempdir()) file.copy(filepaths$base, outdir, recursive = TRUE)


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
  if("character" %in% class(tbl) & length(tbl) != 1) {
    stop("A filepath has to a character vector of length 1")
  }

    if("character" %in% class(tbl)) {

    tbl <- arrow::read_delim_arrow(tbl, ...)

  } else if("data.frame" %in% class(tbl)) {

    tbl <- dplyr::tibble(tbl)

  } else {

    stop("tbl is not a character vector or a data.frame")
  }

  if(!"rowid" %in% colnames(tbl)) tbl$rowid <- 1:nrow(tbl)

  tbl
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------




#' Update CHR/POS/RSID/EffectAllele/OtherAllele for GWAS sumstats using dbSNP
#'
#' @inheritParams tidyGWAS
#' @param data_list a list containing three [dplyr::tibble()]s: main, indel, without_rsid
#' @param filepaths a `list()` of filepaths, as created by [setup_pipeline_paths()]
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#'
#' validate_with_dbsnp(sumstats, build = "NA", dbsnp_path = "/dbsnp155/dbsnp")
#' }
#'
validate_with_dbsnp <- function(data_list, build = c("NA", "37", "38"),filepaths, dbsnp_path, add_missing_build=TRUE) {

  rlang::check_required(dbsnp_path)
  build <- rlang::arg_match(build)


  # -------------------------------------------------------------------------


  cli::cli_h3("7a) Starting with main rows: ")
  if(!is.null(data_list$main)) {
    data_list$main <- repair_dbnsp(data_list$main, dbsnp_path = dbsnp_path, build = build, add_missing_build = add_missing_build)
    filter_func <- make_callback(id = paste0(filepaths$removed_rows, "main_validate_with_dbsnp"))
    data_list$main <- filter_func(data_list$main)

  }

  # handle edge case of length 0 tibble in without RSID
  if(!is.null(data_list$without_rsid)) {
    if(nrow(data_list$without_rsid) == 0) return(dplyr::bind_rows(data_list$main, data_list$indels))
  }

  if(!is.null(data_list$without_rsid)) {
    cli::cli_h3("7b) rows without RSID: ")
    data_list$without_rsid <- repair_dbnsp(data_list$without_rsid, dbsnp_path = dbsnp_path, build = build, add_missing_build = add_missing_build)
    filter_func <- make_callback(id = paste0(filepaths$removed_rows, "without_rsid_validate_with_dbsnp"))
    data_list$without_rsid <- filter_func(data_list$without_rsid)

    # edge cases: -------------------------------------------------------------
    # it is possible that without_rsid subset contains rows that map to the same
    # rsid as already exists in main. Need to handle this

    data_list$main <- dplyr::bind_rows(data_list$main, data_list$without_rsid)
    data_list$without_rsid <- NULL
    before_unique_check <- dplyr::select(data_list$main, rowid)

    # handle duplications
    b38_missing <- dplyr::filter(data_list$main, is.na(CHR)) |>
      dplyr::distinct(CHR_37,POS_37,EffectAllele,OtherAllele, .keep_all = TRUE)

    data_list$main <- dplyr::filter(data_list$main, !is.na(CHR)) |>
      dplyr::distinct(CHR,POS,EffectAllele,OtherAllele, .keep_all = TRUE) |>
      dplyr::bind_rows(b38_missing)


    removed <- dplyr::anti_join(before_unique_check, data_list$main, by = "rowid")

    if(nrow(removed) > 0) {
      cli::cli_alert_info("Found {nrow(removed)} rows which map to a CHR:POS:REF:ALT that another variant maps to. These are removed")
      arrow::write_parquet(removed, paste0(filepaths$removed_rows, "validate_with_dbsnp_duplications.parquet"))
    }


  }

  dplyr::bind_rows(data_list$main, data_list$indels)


}




repair_dbnsp <- function(tbl, dbsnp_path, build, add_missing_build) {
  # existence of chr:pos or rsid decides which columns to repair
  has_rsid <- "RSID" %in% colnames(tbl)
  has_chr_pos <- all(c("CHR", "POS") %in% colnames(tbl))

  if(has_rsid & !has_chr_pos) {
    cli::cli_inform("Using RSID to align with dbSNP")
    tbl <- repair_chr_pos(tbl, dbsnp_path = dbsnp_path, add_missing_build = add_missing_build)

  } else if(has_chr_pos & !has_rsid) {
    cli::cli_inform("Using CHR and POS to align with dbSNP")
    tbl <- repair_rsid(tbl, build = build, dbsnp_path =  dbsnp_path, add_missing_build =add_missing_build)


  } else if(has_chr_pos & has_rsid) {
    cli::cli_inform("Using CHR, POS and RSID to align with dbSNP")
    tbl <- verify_chr_pos_rsid(tbl, build = build, dbsnp_path, add_missing_build = add_missing_build)

  }

  tbl
}



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------




identify_removed_rows <- function(finished, filepaths) {
  start <- arrow::read_parquet(paste0(filepaths$base , "/raw_sumstats.parquet"), select = "rowid")

  removed_rows <- dplyr::anti_join(start, finished, by = "rowid")

  # remove the file with all rowids
  files_in_dir <- list.files(paste0(filepaths$base, "/pipeline_info/"), pattern = "*removed_row*", full.names = TRUE)


  removed_rows_with_flags <-
    files_in_dir |>
    purrr::set_names(base::basename) |>
    purrr::map(arrow::read_parquet) |>
    purrr::map(\(x) dplyr::select(x, rowid)) |>
    purrr::list_rbind(names_to = "reason") |>
    dplyr::select(rowid, reason) |>
    dplyr::mutate(reason = stringr::str_remove(reason, "removed_rows_")) |>
    dplyr::mutate(reason = stringr::str_remove(reason, ".parquet"))

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



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------



#' Create the folder structure for tidyGWAS
#'
#' @inheritParams tidyGWAS
#'
#' @return a list of filepaths
#' @export
#'
#' @examples
#' setup_pipeline_paths("tidyGWAS_first")
setup_pipeline_paths <- function(name) {

  # define workdir
  stopifnot("name for setup_pipeline_paths has to be a character" = is.character(name))
  workdir <- paste(tempdir(), name, name, sep = "/")
  pipeline_info <- paste(workdir,"pipeline_info", sep = "/")
  if(!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
  if(!dir.exists(pipeline_info)) dir.create(pipeline_info, recursive = TRUE)


  list(
    "base" = workdir,
    "logfile" = paste0(workdir, "/tidyGWAS_logfile.txt"),
    "cleaned" = paste(workdir, "tidyGWAS_hivestyle", sep = "/"),

    "updated_rsid"= paste(pipeline_info, "updated_rsid.parquet", sep = "/"),
    "failed_rsid_parse" = paste(pipeline_info, "removed_failed_rsid_parse.parquet", sep = "/"),
    "removed_rows"= paste(pipeline_info, "removed_rows_", sep = "/")
  )





}


write_finished_tidyGWAS <- function(df, output_format, outdir, filepaths) {

  df <- standardize_column_order(df)

  if(output_format == "hivestyle") {

    arrow::write_dataset(dplyr::group_by(df, CHR), filepaths$cleaned)

  } else if(output_format == "parquet") {

    arrow::write_parquet(df, paste(filepaths$base, "cleaned_GRCh38.parquet", sep = "/"))

  } else if(output_format == "csv") {

    outfile <- paste(filepaths$base, "cleaned_GRCh38.csv", sep = "/")
    arrow::write_csv_arrow(df, paste(filepaths$base, "cleaned_GRCh38.csv", sep = "/"))
    system(glue::glue("gzip {outfile}"))

  }


}




# -------------------------------------------------------------------------


#' Create a [dplyr::tibble()] with tidyGWAS column names
#'
#' tidyGWAS functions assumes fixed column names. This function facilitates
#' renaming column into tidyGWAS format.
#'
#' @param tbl a [dplyr::tibble()] or something coercible to one
#' @param CHR chromosome
#' @param POS position
#' @param RSID rsID from dbSNP
#' @param EffectAllele allele corresponding to effect, B
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
      ,"CHR_37", "POS_37", "rowid", "multi_allelic", "indel", "REF" = "ref_allele"))
    )

}


row_removal_report <- function(filepath) {

}

