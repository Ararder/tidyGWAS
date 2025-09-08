library(tidyverse)

# setup_local_dbsnp -------------------------------------------------------
dbsnp_path <- "~/Downloads/dbSNP155_v2/" # on my local macbook
first <- dplyr::select(test_sumstat, CHR, POS, RSID) |>
  dplyr::tibble() |>
  dplyr::mutate(CHR = as.character(CHR))

second <- dplyr::select(pval_as_char_df, CHR, POS, RSID)
second$CHR <- as.character(second$CHR)

data <- dplyr::bind_rows(first, second) |>
  dplyr::distinct() |>
  dplyr::tibble()

grch37_rsid <- map_to_dbsnp(tbl = data, by = "rsid", dbsnp_path = dbsnp_path)
grch37_chr_pos <- map_to_dbsnp(tbl = data, build = "37", by = "chr:pos",dbsnp_path)

final <- dplyr::bind_rows(grch37_chr_pos,grch37_rsid) |>
  dplyr::distinct() |>
  dplyr::mutate(RSID = as.integer(stringr::str_sub(RSID, start = 3)))


arrow::write_dataset(dplyr::group_by(final, CHR), test_path("fixtures/dbSNP"))





# make the dummy refsnp-merged ---------------------------------------------




load(test_path("data/test_sumstat.rds"))
load(test_path("data/b38_t1d_chr_pos_rsid_pvalue_as_character.rds"))

temp <- dplyr::bind_rows(
  dplyr::select(test_file, RSID),
  dplyr::select(pval_as_char_df, RSID)
) |>
  dplyr::mutate(rowid = 1:127715) |>
  flag_rsid_history()

dir.create(test_path("fixtures/dbSNP155/refsnp-merged/"))
dplyr::filter(temp, !is.na(old_RSID)) |>
  dplyr::select(old_RSID, RSID) |>
  arrow::write_parquet(test_path("fixtures/dbSNP155/refsnp-merged/part-0.parquet"))




# -------------------------------------------------------------------------
# make EAF file
df <- arrow::read_parquet("~/Downloads/HRC_eur_0.001.parquet")
all <- dplyr::select(pval_as_char_df, RSID) |>
  dplyr::bind_rows(
    dplyr::select(test_sumstat, RSID)
  ) |>
  dplyr::distinct()

sub <- df |>
  dplyr::filter(RSID %in% all$RSID)

arrow::write_parquet(sub, test_path("fixtures/HRC_eur_0.001.parquet"))



#  make EAF V2 ------------------------------------------------------------
df <- arrow::open_dataset(test_path("fixtures/dbSNP155/v155")) |> dplyr::collect()
all <- arrow::open_dataset("~/Downloads/dbSNP155/EAF_REF_1KG/")

colnames(df)
all |>
  dplyr::mutate(POS = as.integer(POS)) |>
  dplyr::semi_join(df, by = c("CHR" ="CHR", "POS" = "POS_38")) |>
    dplyr::group_by(ancestry) |>
    arrow::write_dataset(test_path("fixtures/dbSNP155/EAF_REF_1KG/"))


arrow::open_dataset(test_path("fixtures/dbSNP155/EAF_REF_1KG/")) |> dplyr::collect()

# make ancestry file
rs <- dplyr::bind_rows(tbl, dplyr::select(pval_as_char_df, -CHR, -P)) |>
  dplyr::distinct(RSID)
df <- arrow::read_parquet("~/Downloads/dbSNP155/ancestry_data.parquet")
out <- dplyr::semi_join(df, rs)

arrow::write_parquet(out, test_path("fixtures/dbSNP155/ancestry_data.parquet"))






# make indels data
dbsnp_path <- "~/Downloads/dbSNP155/"
kk <- bind_rows(select(pval_as_char_df, CHR, POS, EffectAllele, OtherAllele)) |> mutate(CHR = as.character(CHR))
matches <- map_indels_dbsnp(kk, dbsnp_path, by = "chr:pos")

kk2 <- bind_rows(select(test_sumstat, CHR, POS, EffectAllele, OtherAllele)) |> mutate(CHR = as.character(CHR))
matches2 <- map_indels_dbsnp(kk2, dbsnp_path, by = "chr:pos")

all_variants <- bind_rows(matches, matches2)



arrow::open_dataset(paste0(dbsnp_path, "dbSNP157_indels", "/37"))
bind_rows()
arrow::write_dataset()

test
all_variants |>
  select(CHR, POS = POS_38, RSID, ALT = EffectAllele, REF = REF_38) |>
  arrow::write_dataset(test_path("fixtures/dbSNP155/dbSNP157_indels/38/"))


all_variants |>
  select(CHR, POS = POS_37, RSID, ALT = EffectAllele, REF = REF_37) |>
  arrow::write_dataset(test_path("fixtures/dbSNP155/dbSNP157_indels/37/"))




