load(test_path("data/test_sumstat.rds"))
test_sumstat <- dplyr::tibble(test_sumstat)
test_sumstat$CHR <- as.character(test_sumstat$CHR)
load(test_path("data/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))
pval_as_char_df$rowid <- 1:nrow(pval_as_char_df)
dbsnp_path <- test_path("fixtures/dbSNP155")



# default args
logfile <- FALSE
name <- "testing"
convert_p <- 0
verbose <- FALSE
overwrite = FALSE
build <- "NA"
tbl <- test_sumstat
output_dir = tempfile()
indel_strategy <- "keep"
default_build <- "37"
allow_duplications <- FALSE
repair_cols <- TRUE


setup_test_data_meta_analysis <- function(EAF = NULL, INFO = NULL, N = NULL) {
  tmp <- test_sumstat |> dplyr::mutate(N = CaseN + ControlN, REF_37 = EffectAllele)
  meta1 <- tmp[1,]

  ten_sumstats <- dplyr::tibble(
    CHR = meta1$CHR,
    POS = meta1$POS,
    RSID = meta1$RSID,
    EffectAllele = meta1$EffectAllele,
    OtherAllele = meta1$OtherAllele,
    B = meta1$B,
    EAF = EAF,
    SE = meta1$SE,
    REF_37 = EffectAllele,
    INFO = INFO,
    N = N
  )

  tmpdir <- tempfile()

  make_fake_sumstats <- function(tbl) {
    tmpdir <- tempdir()
    outdir <- fs::dir_create(fs::path(tmpdir, "fake_sumstats"))
    for(i in 1:nrow(tbl)) {
      out <- paste0(outdir, "/trait", i)
      arrow::write_dataset(dplyr::group_by(tbl |> dplyr::slice(i), CHR), out)
    }
    outdir
  }
  outdir <- make_fake_sumstats(ten_sumstats)

  m <- fs::dir_ls(outdir) |>
    purrr::map(arrow::open_dataset) |>
    purrr::reduce(c)

  m
}
