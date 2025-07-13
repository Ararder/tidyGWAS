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

  base_variants <- dplyr::slice_sample(test_sumstat, n = 5, by = CHR) |>
    dplyr::select(
      CHR, POS, RSID, EffectAllele, OtherAllele
    )

  var <- base_variants
  permute_sumstats <- function(
    var, EAF = NULL, INFO = NULL, N = NULL,
    EffectiveN=NULL, CaseN=NULL, ControlN=NULL,
    B = rnorm(110, 0, 0.1),
    SE = rnorm(110, 0.05, 0.01)

    ) {


    dplyr::mutate(
      var, EAF = EAF, INFO = INFO, N = N, EffectiveN = EffectiveN,
      CaseN = CaseN, ControlN = ControlN, B=B, SE=SE)
  }

  permute_sumstats(base_variants, EAF = runif(110))




  ten_sumstats <- dplyr::tibble(
    CHR = meta1$CHR,
    POS_37 = meta1$POS,
    POS_38 = meta1$POS,
    RSID = meta1$RSID,
    EffectAllele = meta1$EffectAllele,
    OtherAllele = meta1$OtherAllele,
    B = meta1$B,
    EAF = EAF,
    SE = meta1$SE,
    REF_37 = EffectAllele,
    REF_38 = EffectAllele,
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



setup_test_data_meta_analysis <- function(
    test_sumstat,
    n_variants_per_chr   = 5,
    n_traits             = 10,
    frac_missing_N       = 0,
    frac_missing_EAF     = 0,
    frac_missing_INFO    = 0,
    drop_N_for_traits    = FALSE,
    drop_EAF_for_traits  = FALSE,
    drop_INFO_for_traits = FALSE,
    seed = NULL,
    outdir = tempfile()
) {



  # Recycle the *drop* logicals to length `n_traits` ---------------------------
  recycle_logical <- function(x) {
    if (length(x) == 1L) return(rep.int(x, n_traits))
    if (length(x) != n_traits) {
      stop("Logical vector must be length 1 or `n_traits`.")
    }
    as.logical(x)
  }
  drop_N_for_traits    <- recycle_logical(drop_N_for_traits)
  drop_EAF_for_traits  <- recycle_logical(drop_EAF_for_traits)
  drop_INFO_for_traits <- recycle_logical(drop_INFO_for_traits)

  # ---- 2.  Sample the base variant panel ------------------------------------
  base_variants <- test_sumstat |>
    dplyr::group_by(CHR)                             |>
    dplyr::slice_head(n = n_variants_per_chr)      |>
    dplyr::ungroup() |>
    dplyr::mutate(N = CaseN + ControlN) |>
    dplyr::select(CHR, POS, RSID, EffectAllele, OtherAllele) |>
    dplyr::mutate(POS_37 = POS, POS_38=POS)

  n_variants <- nrow(base_variants)

  # Helper to apply per-row missingness ---------------------------------------
  introduce_na <- function(vec, fraction) {
    if (fraction <= 0) return(vec)
    n_na <- ceiling(fraction * length(vec))
    vec[sample(length(vec), n_na)] <- NA
    vec
  }

  # ---- 3.  Generate trait-specific tables -----------------------------------
  trait_tbls <- purrr::map(seq_len(n_traits), function(i) {
    tbl <- base_variants |>
      dplyr::mutate(
        Trait       = paste0("trait", i),
        B           = stats::rnorm(n_variants, 0, 0.10),
        SE          = stats::rnorm(n_variants, 0.05, 0.01),
        N           = sample(5e3:3e4, n_variants, replace = TRUE),
        EAF         = stats::runif(n_variants, 0.01, 0.99),
        INFO        = stats::runif(n_variants, 0.3,  1.00),
        CaseN       = sample(2e3:1e4, n_variants, replace = TRUE),
        ControlN    = sample(2e3:1e4, n_variants, replace = TRUE),
        EffectiveN  = CaseN + ControlN,
        REF_37      = EffectAllele,
        REF_38      = EffectAllele
      )

    # Column-level missingness -------------------------------------------------
    if (!drop_N_for_traits[i])
      tbl$N   <- introduce_na(tbl$N,   frac_missing_N) else tbl$N   <- NULL
      if (!drop_EAF_for_traits[i])
        tbl$EAF <- introduce_na(tbl$EAF, frac_missing_EAF) else tbl$EAF <- NULL
        if (!drop_INFO_for_traits[i])
          tbl$INFO <- introduce_na(tbl$INFO, frac_missing_INFO) else tbl$INFO <- NULL

          tbl
  })


  purrr::walk(seq_along(trait_tbls), function(i) {
    arrow::write_dataset(
      dplyr::group_by(trait_tbls[[i]], CHR),
      fs::path(outdir, paste0("trait", i))
    )
  })

  # ---- 5.  Return a lazily concatenated Arrow Dataset -----------------------
  arrow::open_dataset(outdir)
}

