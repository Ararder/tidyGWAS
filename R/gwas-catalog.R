
from_gwas_catalog <- function(study_id, quiet=FALSE) {


  rlang::is_scalar_character(study_id) || cli::cli_abort("study_id must be a single character string")
  is_sumstats_available(study_id) || cli::cli_abort("Full summary statistics for {study_id} must be available in GWAS catalog")
  workdir <- fs::dir_create(fs::path(tempdir(), study_id))


  all_urls <- get_urls(study_id)
  base_gwas <- all_urls[stringr::str_detect(all_urls, "tsv.gz$")]
  base_meta <- all_urls[stringr::str_detect(all_urls, "meta.yaml$")]
  meta_file <- fs::path(workdir, fs::path_file(base_meta))
  gwas_file <- fs::path(workdir, fs::path_file(base_gwas))


  # download meta file and print it?
  curl::curl_download(url = base_meta, destfile = meta_file)
  meta <- yaml::read_yaml(meta_file)


  cli::cli_alert_info("Downloading GWAS summary statistics → {.path {gwas_file}}")
  curl::curl_download(url = base_gwas,destfile = gwas_file, quiet=quiet)

  gwas_file
}


is_sumstats_available <- function(study_id, server = "https://www.ebi.ac.uk/gwas/rest/api") {


  resp <- httr::GET(paste0(server,"/studies/", study_id), httr::accept("application/json"))
  httr::stop_for_status(resp)
  jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))$fullPvalueSet

}

get_urls <- function(study_id) {
  root <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

  num <- as.integer(sub("GCST", "", study_id))
  lo  <- ( (num - 1) %/% 1000 ) * 1000 + 1
  hi  <- lo + 999

  range_folder <- sprintf("GCST%08d-GCST%08d", lo, hi)
  subdir <- file.path(root, range_folder, study_id)
  harmon <- file.path(subdir, "harmonised")

  scrape_dir(subdir)
}


scrape_dir <- function(url) {
  page  <- xml2::read_html(url)
  hrefs <- rvest::html_attr(rvest::html_nodes(page, "a"), "href")
  hrefs <- hrefs[!hrefs %in% c("../", "", NA)]
  file.path(url, hrefs)
}





