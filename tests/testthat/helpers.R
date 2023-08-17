randomizer <- function(end=10) {
  sample(c(1:10),nrow(tbl), replace=TRUE) == 1
}
mock_dbsnp <- function(){
  # set .env to parent.frame() so the mockings persist in the parent function
  load(test_path("fixtures/b38.rds"))
  load(test_path("fixtures/b37.rds"))
  load(test_path("fixtures/rs_merge_arch.rds"))

  local_mocked_bindings(
    snpsById = function(x, ids, ...) {
      if(x == 38) {
        dplyr::filter(b38, RSID %in% ids) |>
          dplyr::rename(seqnames = CHR, pos = POS, RefSNP_id = RSID)
      } else if(x == 37) {
        dplyr::filter(b37, RSID %in% ids) |>
          dplyr::rename(seqnames = CHR, pos = POS, RefSNP_id = RSID)
      }

    },
    .env = parent.frame(),
    .package = "BSgenome"
    )

  local_mocked_bindings(
    snpsByOverlaps = function(x, ranges, ...) {
      if(x == 38) {
        dplyr::semi_join(b38, ranges, by = c("CHR" = "seqnames", "POS" = "pos")) |>
          dplyr::rename(seqnames = CHR, pos = POS, RefSNP_id = RSID)
      } else if(x == 37) {
        dplyr::semi_join(b37, ranges, by = c("CHR" = "seqnames", "POS" = "pos")) |>
          dplyr::rename(seqnames = CHR, pos = POS, RefSNP_id = RSID)
      }

    },
    .env = parent.frame(),
    .package = "BSgenome"
  )


  local_mocked_bindings(
    get_bsgenome = function() {
      list("snps_37" = 37, "snps_38" = 38, "genome_37" = "genome", "genome_38" = "genome")
    },
    .package = "tidyGWAS",
    .env = parent.frame(),
  )

  local_mocked_bindings(
    get_ref_data = function() {
      rs_merge_arch
    },
    .env = parent.frame(),
  )
}

mock_arrow <- function() {
  withr::local_envvar(
    "rs_merge_arch" = test_path("fixtures/RsMergeArch.parquet"),
    .local_envir = parent.frame()
    )

  withr::local_envvar(
    "grch37" = test_path("fixtures/dbSNP/GRCh37"),
    .local_envir = parent.frame()
  )
  withr::local_envvar(
    "grch38" = test_path("fixtures/dbSNP/GRCh38"),
    .local_envir = parent.frame()
  )

}








