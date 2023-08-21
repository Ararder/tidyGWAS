utils::globalVariables(c(
  "chr_pos", "n", "n_chr_pos", "dup_rsid",
  ".data", "invalid_rsid",  "new_rsid",
  "has_rsid", "head", "CHR_37", "POS_37",
   "reason", "filter_callback", "no_dbsnp_entry", "logfile"))

snp_cols <- c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele", "rowid")
info_cols <- c("INFO", "N", "CaseN", "ControlN", "EAF")
stats_cols <- c("B", "Z", "OR", "P", "SE", info_cols)
valid_column_names <- c(snp_cols, stats_cols, info_cols)




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
#' @param ... pass additional arguments to [arrow::read_delim_arrow()] which is the function
#' that will be called if tbl is a filepath
#' possible arguments are `study_n` to set N, and `build` to set genome build.
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
#' @param study_n Use to set N in tbl
#' @param keep_indels Should indels be kept?
#' @param repair_cols Should any missing columns be repaired?
#' @param logfile Should messages be redirected to a logfile?
#' @param log_on_err Optional. Can pass a filepath to copy the logfile to when the function exists.
#' @param verbose Explain filters in detail?
#'
#' @return a tibble or NULL, depending on outdir
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
    convert_p=2.225074e-308,
    name = stringr::str_replace_all(date(), pattern = c(" "="_", ":"="_")),
    keep_indels = TRUE,
    repair_cols = TRUE,
    logfile = FALSE,
    log_on_err="tidyGWAS_logfile.txt",
    verbose = FALSE
    ) {


  # parse arguments ---------------------------------------------------------
  tbl <- parse_tbl(tbl, ...)

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
  start_time <- Sys.time()
  rows_start <- nrow(tbl)

  # welcome message ----------------------------------------------------------

  cli::cli_h1("Running {.pkg tidyGWAS {packageVersion('tidyGWAS')}}")
  cli::cli_inform("Starting at {start_time}, with {rows_start} in input data.frame")
  cli::cli_alert_info("Saving all files during execution to {.file {filepaths$base}},")
  cli::cli_alert_info("After execution, files will be copied to {.file {paste(outdir,name, sep = '/')}}")
  cli::cli_alert_info("Detected {rows_start} rows in input summary statistics")


  tbl <- select_correct_columns(tbl, study_n)

  cli::cli_h2("1) Scanning for rows with NA")
  tbl <- remove_rows_with_na(tbl, filepaths)

  cli::cli_h2("2) Updating merged RSIDs")
  if("RSID" %in% colnames(tbl) & !missing(dbsnp_path)) tbl <- update_rsid(tbl, dbsnp_path = dbsnp_path, filepaths = filepaths)

  cli::cli_h2("3) Scanning for rows with duplications")
  tbl <- remove_duplicates(tbl, filepaths = filepaths)

  # setup list
  data_list <-  list("main" = NULL, "indels" = NULL, "without_rsid" = NULL)
  cli::cli_h2("4) Scanning for indels")
  tmp <- detect_indels(tbl, keep_indels)
  data_list$main <- tmp$main
  if(!is.null(tmp$indels)) data_list$indels <- tmp$indels



  if("RSID" %in% colnames(tbl)) {
    cli::cli_h2("5) Scanning for invalid RSID ")
    unpack <- validate_rsid(data_list$main, outpath = paste(filepaths$removed_rows, "invalid_chr_pos_rsid.parquet"))
    data_list$main <- unpack$main
    if(!is.null(unpack$without_rsid)) data_list$without_rsid <- unpack$without_rsid
  }

  # there is three possible subsets of sumstats now: main, indels, and without CHR:POS
  # apply column validations on each of them
  filter_funcs <-  purrr::map(paste0(filepaths$removed_rows, c("main", "indels", "without_rsid"), "_col_validation"), make_callback)
  cols_to_not_validate <- list("", c("EffectAllele","OtherAllele"), "")
  id <- list("main rows", "indel rows", "rows without RSID")

  data_list <- list(tbl = data_list, remove_cols = cols_to_not_validate, filter_func = filter_funcs, verbose = list(verbose), convert_p = list(convert_p), id = id) |>
    purrr::pmap(validate_sumstat)



  # fix CHR/POS/RSID --------------------------------------------------------


  if(!missing(dbsnp_path)) {

    filter_funcs <-  purrr::map(paste0(filepaths$removed_rows,  c("main", "without_rsid"), "_incompat_alleles_or_not_in_dbsnp"), make_callback)
    cli::cli_h2("Validating CHR/POS/RSID dbSNP for main rows")
    main <- validate_with_dbsnp(data_list$main, build = build, dbsnp_path, filter_func = filter_funcs[[1]])

    if(!is.null(data_list$without_rsid)) {

      cli::cli_h2("Validating CHR/POS/RSID using dbSNP for rows without RSID")
      data_list$without_rsid <- validate_with_dbsnp(data_list$without_rsid, build = build, dbsnp_path, filter_func = filter_funcs[[2]])

      main <- dplyr::bind_rows(main,data_list$without_rsid)

      # handle duplications
      b38_missing <- dplyr::filter(main, is.na(CHR)) |>
        dplyr::distinct(CHR_37,POS_37,EffectAllele,OtherAllele, .keep_all = TRUE)

      main <- dplyr::filter(main, !is.na(CHR)) |>
        dplyr::distinct(CHR,POS,EffectAllele,OtherAllele, .keep_all = TRUE) |>
        dplyr::bind_rows(b38_missing)



    }
  } else {
    main <- data_list$main
  }


  # merge together and repair missing statistics columns ---------------------

  if(keep_indels) main <- dplyr::bind_rows(main, data_list$indels)
  if(repair_cols) main <- repair_stats(main)

  # flag indels and multi-allelics
  main <- flag_duplicates(main, "rsid") |>
    flag_indels() |>
    dplyr::rename(multi_allelic = dup_rsid)
  identify_removed_rows(dplyr::select(main,rowid), filepaths)


  # end of pipeline  --------------------------------------------------------

  fmt <- prettyunits::pretty_dt(Sys.time() - start_time)
  cli::cli_h1("Finished tidyGWAS")
  cli::cli_alert_info("A total of {rows_start - nrow(main)} rows were removed")
  cli::cli_alert_info("Total running time: {fmt}")
  write_finished_tidyGWAS(df = main, output_format = output_format, outdir = outdir, filepaths = filepaths)

  }


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



#' Update CHR/POS/RSID/EffectAllele/OtherAllele for GWAS sumstats using dbSNP
#'
#' @param tbl a [dplyr::tibble()] with column names in [tidyGWAS_columns()] format.
#' @param build Genome build of sumstats. if `'NA'` [infer_build()] will be used
#' to automatically detect build.
#' @param dbsnp_path filepath to dbSNP files in .parquet format.
#' @param filter_func A function that is applied to the tibble before exiting. This is used
#' to write to disk any rows that removed.
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#'
#' validate_with_dbsnp(sumstats, build = "NA", dbsnp_path = "/dbsnp155/dbsnp")
#' }
#'
validate_with_dbsnp <- function(tbl, build = c("NA", "37", "38"), dbsnp_path, filter_func) {
  build <- rlang::arg_match(build)
  if(nrow(tbl) == 0) return(tbl)
  cli::cli_li("Remove rows where REF/ALT in dbSNP is not compatible with EffectAllele / OtherAllele")
  cli::cli_li("Remove rows where RSID/CHR:POS does not match a entry in dbSNP v.155")

  has_rsid <- "RSID" %in% colnames(tbl)
  has_chr_pos <- all(c("CHR", "POS") %in% colnames(tbl))

  # existence of chr:pos or rsid decides which columns to repair
  if(has_rsid & !has_chr_pos) {

    cli::cli_li("Repairing chromosome and position")
    main_df <- repair_chr_pos(tbl, dbsnp_path = dbsnp_path)

  } else if(has_chr_pos & !has_rsid) {

    cli::cli_li("Repairing RSID")
    main_df <- repair_rsid(tbl, build = build, dbsnp_path =  dbsnp_path)


  } else if(has_chr_pos & has_rsid) {

    cli::cli_li("Checking that CHR:POS and RSID match. RSID will be updated accordingly to dbSNP")
    main_df <- verify_chr_pos_rsid(tbl, build = build, dbsnp_path)

  }


  if(!missing(filter_func)) main_df <- filter_func(main_df)

  main_df



}







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


  if(outdir != tempdir()) {

    file.copy(filepaths$base, outdir, recursive = TRUE)

  }

  df
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
      ,"CHR_37", "POS_37", "rowid", "multi_allelic", "indel"))
    )

}
