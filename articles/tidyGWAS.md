# tidyGWAS

## Installation

tidyGWAS is not yet on CRAN. To install, you need to use
`devtools::install_github()` or `remotes::install_github()`

``` r
devtools::install_github("ararder/tidyGWAS")
remotes::install_github("ararder/tidyGWAS")
```

Secondly, you will need to download the reference data. tidyGWAS uses a
slightly edited version of dbSNP v155, that is available
[here](https://zenodo.org/records/16639374).

``` bash
wget https://zenodo.org/records/16639374/files/dbSNP155.tar
tar -xvf dbSNP155.tar
```

### Using an apptainer or docker container instead

You can skip the installation step by using a docker container with
docker or apptainer.

``` bash
# using apptainer
apptainer pull docker://arvhar/tidygwas:latest
apptainer shell ~/tidygwas_latest.sif

# using docker 
docker run arvhar/tidygwas:latest
```

## NEWS

With the release of **tidyGWAS 1.0**, the package now supports automatic
column-name guessing as well as automatic detection of file delimiter,
using
[`data.table::fread()`](https://rdatatable.gitlab.io/data.table/reference/fread.html).
This means that the `delim` argument has been superseded and it’s often
no longer nessecary to pass the names of the input column names with the
`column_names` argument. But remember to double check the logfile that
the column names were parsed correctly.

## Quick start

This is what a typical call to
[`tidyGWAS()`](https://ararder.github.io/tidyGWAS/reference/tidyGWAS.md)
will look like. A detailed explanation of the different arguments
follows.

``` r
# we setup a directory where we will store all summary statistics cleaned by tidyGWAS
gwas_folder <- tempfile()
# provide the filepath to the name of a directory that tidyGWAS will create.
outdir <- paste0(gwas_folder, "gwasName")

cleaned <- tidyGWAS(
  tbl = "/filepath/to/sumstats", 
  dbsnp_path = dsnp_path,
  # Number of samples are missing, so manually impute
  CaseN = 54000,
  ControlN = 73000,
  logfile=TRUE,
  output_dir = outdir
  )
```

## Tutorial

tidyGWAS now supports directly downloading files from GWAS catalog or
providing URLs.

The first argument to tidyGWAS is `tbl`, and should be one off

1.  A [GWAS catalog](https://www.ebi.ac.uk/gwas/) study id. For example,
    “[GCST90101808](https://www.ebi.ac.uk/gwas/studies/GCST90101808)”.
    The full summary statistics needs to be available.
2.  A webpage URL to a downloadable file. For example, the FinnGEN R[12
    depression](https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_F5_DEPRESSIO.gz)
    release
3.  A local filepath. For example “~/summary_stats/raw/trait_x.tsv.gz”
4.  An in-memory data.frame.

``` r
library(tidyGWAS)
# we use the dummy version of dbSNP that comes with the package
dbsnp_path <- system.file("extdata/dbSNP155", package = "tidyGWAS")

# a dummy sumstats with 100 000 rows
gwas <- system.file("extdata/sumstats.tsv.gz", package = "tidyGWAS")
# store the results in a temporary directory
out <- tempfile()

# using a GWAS catalog study ID
tidyGWAS(
  tbl = "GCST90101808",
  dbsnp_path = dbsnp_path,
  output_dir = "finnGEN_MDD_R12"
)

# using a URL
tidyGWAS(
  tbl = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_F5_DEPRESSIO.gz",
  dbsnp_path = dbsnp_path,
  output_dir = "finnGEN_MDD_R12"
)

# using a local flepath
gwas <- system.file("extdata/sumstats.tsv.gz", package = "tidyGWAS")
tidyGWAS(
  tbl = gwas,
  dbsnp_path = dbsnp_path,
  output_dir = "test_sumstat"
)

# or an in-memory data.frame
df = read.table(gwas)
tidyGWAS(
  tbl = df,
  dbsnp_path = dbsnp_path,
  output_dir = "test_sumstat"
)
```

### An example call

Here’s what the output of tidyGWAS looks like, using the test file
provided with tidyGWAS

``` r
library(tidyGWAS)
tidyGWAS(
  tbl = system.file("extdata/sumstats.tsv.gz", package = "tidyGWAS"),
  dbsnp_path = system.file("extdata/dbSNP155", package = "tidyGWAS"),
  output_dir = tempfile()
)
#> Parsing input summary statistics...
#> 
#> 
#> ── Running tidyGWAS 1.0.0 ──────────────────────────────────────────────────────
#> 
#> Starting at 2025-11-29 11:34:14.270319
#> with 100000 rows in input data.frame
#> ℹ Saving output in folder: /tmp/RtmpRYimIr/file1b62d057865
#> 
#> 
#> 
#> ── The detected format is "tidyGWAS" 
#> 
#> Was able to map 12 out of 12 to the "tidyGWAS" format
#> ✔ CHR -> CHR
#> 
#> ✔ POS -> POS
#> 
#> ✔ RSID -> RSID
#> 
#> ✔ EffectAllele -> EffectAllele
#> 
#> ✔ OtherAllele -> OtherAllele
#> 
#> ✔ B -> B
#> 
#> ✔ SE -> SE
#> 
#> ✔ EAF -> EAF
#> 
#> ✔ P -> P
#> 
#> ✔ CaseN -> CaseN
#> 
#> ✔ ControlN -> ControlN
#> 
#> ✔ INFO -> INFO
#> 
#> 
#> 
#> ── Checking that columns follow tidyGWAS format 
#> 
#> ✔ The following columns are used for further steps: CHR, POS, RSID, EffectAllele, OtherAllele, B, SE, EAF, INFO, P, CaseN, ControlN, and rowid
#> 
#> 
#> 
#> ── Checking for columns with all NA 
#> 
#> ✔ Found no columns with all NA
#> 
#> ℹ Found CaseN and ControlN, and no effective N:
#> Calculating EffectiveN by `EffectiveN = 4 / (1 / ControlN + 1 / CaseN)`
#> 
#> ℹ Found CaseN and ControlN, Calculating N by `N = ControlN + CaseN`
#> 
#> 
#> 
#> ── 1) Scanning for rows with NA in critical columns ──
#> 
#> 
#> 
#> ✔ No rows contained missing values in CHR, POS, EffectAllele, and OtherAllele
#> 
#> 
#> 
#> ── 2) Scanning for rows with duplications ──
#> 
#> 
#> 
#> ℹ Looking for duplications with columns: CHR, POS, EffectAllele, and OtherAllele
#> 
#> ✔ Found no duplications
#> 
#> 
#> 
#> ── 3) Scanning for indels ──
#> 
#> 
#> 
#> 1. EffectAllele or OtherAllele, character length > 1: A vs AA
#> 
#> 2. EffectAllele or OtherAllele coded as 'D', 'I', or 'R'
#> 
#> ✔ Detected 0 rows as indels
#> 
#> 
#> 
#> ── 4) Validating columns 
#> 
#> ℹ The median value of B is -0.0012, which seems reasonable
#> 
#> 
#> 
#> ── All rows passed validation 
#> 
#> 
#> 
#> ── 5) Adding RSID based on CHR:POS. Adding dbSNP based QC flags 
#> 
#> ℹ Inferring build by matching 10000 rows to GRCh37 and GRCh38
#> 
#> 99 snps matched GRCh38, 9998 for GRCh37, inferring build to be 37
#> ℹ 0 rows had EffectAllele and OtherAllele not matching dbSNP REF/ALT alleles,but which matched after strand-flipping
#> 
#> ! Removed 21 rows with no dbSNP entry or with incompat alleles
#> 
#> /tmp/RtmpRYimIr/file1b62d057865/pipeline_info/removed_nodbsnp.parquet
#> 
#> 
#> ── 6) Repairing missings statistics columns if possible ──
#> 
#> 
#> 
#> ℹ Z missing: Calculating Z using the formula:  Z = B / SE
#> 
#> 
#> 
#> ── Finished repair_stats:  
#> 
#> ℹ Added 1 new columns: Z
#> 
#> 
#> 
#> ── Listing final breakdown of removed rows:  
#> 
#> nodbsnp: 21
#> 
#> 
#> 
#> ── Finished tidyGWAS ───────────────────────────────────────────────────────────
#> 
#> ℹ A total of 21 rows were removed
#> 
#> ℹ Total running time: 3.8s
#> 
#> Saving metadata from analysis to /tmp/RtmpRYimIr/file1b62d057865/metadata.yaml
#> # A tibble: 99,979 × 21
#>    CHR   POS_37 EffectAllele OtherAllele rowid     B        P     SE  INFO CaseN
#>    <chr>  <int> <chr>        <chr>       <int> <dbl>    <dbl>  <dbl> <dbl> <int>
#>  1 6     2.75e7 C            T           77548 0.212 4.83e-39 0.0162 0.99  53386
#>  2 6     2.81e7 T            C           36396 0.207 4.67e-38 0.016  1     53386
#>  3 6     2.75e7 C            T           77992 0.206 7.92e-38 0.016  0.993 53386
#>  4 6     2.78e7 T            C           67347 0.171 1.45e-32 0.0144 1     53386
#>  5 6     2.63e7 T            C            7577 0.191 8.94e-30 0.0168 0.996 53386
#>  6 6     3.08e7 G            A           26615 0.157 9.77e-28 0.0144 0.954 53386
#>  7 6     2.83e7 C            T           91092 0.116 4.33e-27 0.0108 1.01  53386
#>  8 6     2.64e7 C            T           84820 0.151 9.52e-27 0.0141 1     53386
#>  9 6     3.11e7 C            T           17683 0.143 1.54e-25 0.0137 1     53386
#> 10 6     3.03e7 A            G           71635 0.154 2.58e-25 0.0148 0.992 49683
#> # ℹ 99,969 more rows
#> # ℹ 11 more variables: ControlN <int>, EAF <dbl>, EffectiveN <int>, N <int>,
#> #   POS_38 <int>, RSID <chr>, REF_37 <chr>, REF_38 <chr>, indel <lgl>,
#> #   multi_allelic <lgl>, Z <dbl>
```

### Manually setting column names

tidyGWAS will attempt to guess the column names of the input file.
Always take a look at the logfile if it’s the first time you are using
that particular format. Sometimes the automatic parsing fails, like in
the example below

``` r
testdf <- arrow::read_tsv_arrow(system.file("extdata/sumstats.tsv.gz", package = "tidyGWAS"))

# nonsense names
names(testdf) <- c("CHROM_X", "BP_NEW", "RsID", "A1", "A2", "B","SE", "EAF","INFO", "P", "CaseN", "ControlN", "rowid")

# Wit
tryCatch(
  tidyGWAS(
    tbl = testdf,
    dbsnp_path = system.file("extdata/dbSNP155", package = "tidyGWAS"),
    output_dir = tempfile()
  ),
  error = function(e) {
    message("tidyGWAS failed: ", e$message)
    NULL
  }
)
#> Parsing input summary statistics...
#> 
#> 
#> ── Running tidyGWAS 1.0.0 ──────────────────────────────────────────────────────
#> 
#> Starting at 2025-11-29 11:34:18.247211
#> with 100000 rows in input data.frame
#> ℹ Saving output in folder: /tmp/RtmpRYimIr/file1b62765d05d
#> 
#> 
#> 
#> ── The detected format is "tidyGWAS" 
#> 
#> Was able to map 7 out of 12 to the "tidyGWAS" format
#> ✔ B -> B
#> 
#> ✔ SE -> SE
#> 
#> ✔ EAF -> EAF
#> 
#> ✔ P -> P
#> 
#> ✔ CaseN -> CaseN
#> 
#> ✔ ControlN -> ControlN
#> 
#> ✔ INFO -> INFO
#> 
#> ! Failed to map: "CHROM_X", "BP_NEW", "RsID", "A1", and "A2"
#> 
#> ℹ Attempting to map: "CHROM_X", "BP_NEW", "RsID", "A1", and "A2" to a dictionary of column names
#> 
#> ✔ EffectAllele -> A1
#> 
#> ✔ OtherAllele -> A2
#> 
#> ! extra column(s) "CHROM_X", "BP_NEW", and "RsID" were not mapped to a tidyGWAS column
#> 
#> 
#> 
#> ── Checking that columns follow tidyGWAS format 
#> 
#> ✔ The following columns are used for further steps: EffectAllele, OtherAllele, B, SE, EAF, INFO, P, CaseN, ControlN, and rowid
#> 
#> ✖ Removed columns:  CHROM_X, BP_NEW, and RsID
#> 
#> 
#> 
#> ── Checking for columns with all NA 
#> 
#> ✔ Found no columns with all NA
#> 
#> tidyGWAS failed: Either CHR and POS or RSID are required columns
#> NULL
```

If you use column_names, those columns will be renamed, following by an
automatic attempt at parsing the remaining column names.

``` r
testdf <- arrow::read_tsv_arrow(system.file("extdata/sumstats.tsv.gz", package = "tidyGWAS"))

# nonsense names
names(testdf) <- c("CHROM_X", "BP_NEW", "RsID", "A1", "A2", "B","SE", "EAF","INFO", "P", "CaseN", "ControlN", "rowid")
out <- tempfile()
tidyGWAS(
  tbl = testdf,
  dbsnp_path = system.file("extdata/dbSNP155", package = "tidyGWAS"),
  output_dir = out,
  column_names = list(
    CHR = "CHROM_X",
    POS = "BP_NEW",
    RSID = "RsID"
  )
)
#> Parsing input summary statistics...
#> 
#> 
#> ── Running tidyGWAS 1.0.0 ──────────────────────────────────────────────────────
#> 
#> Starting at 2025-11-29 11:34:18.591913
#> with 100000 rows in input data.frame
#> ℹ Saving output in folder: /tmp/RtmpRYimIr/file1b6234d4c157
#> 
#> 
#> 
#> ── The detected format is "tidyGWAS" 
#> 
#> Was able to map 10 out of 12 to the "tidyGWAS" format
#> ✔ CHR -> CHR
#> 
#> ✔ POS -> POS
#> 
#> ✔ RSID -> RSID
#> 
#> ✔ B -> B
#> 
#> ✔ SE -> SE
#> 
#> ✔ EAF -> EAF
#> 
#> ✔ P -> P
#> 
#> ✔ CaseN -> CaseN
#> 
#> ✔ ControlN -> ControlN
#> 
#> ✔ INFO -> INFO
#> 
#> ! Failed to map: "A1" and "A2"
#> 
#> ℹ Attempting to map: "A1" and "A2" to a dictionary of column names
#> 
#> ✔ EffectAllele -> A1
#> 
#> ✔ OtherAllele -> A2
#> 
#> 
#> 
#> ── Checking that columns follow tidyGWAS format 
#> 
#> ✔ The following columns are used for further steps: CHR, POS, RSID, EffectAllele, OtherAllele, B, SE, EAF, INFO, P, CaseN, ControlN, and rowid
#> 
#> 
#> 
#> ── Checking for columns with all NA 
#> 
#> ✔ Found no columns with all NA
#> 
#> ℹ Found CaseN and ControlN, and no effective N:
#> Calculating EffectiveN by `EffectiveN = 4 / (1 / ControlN + 1 / CaseN)`
#> 
#> ℹ Found CaseN and ControlN, Calculating N by `N = ControlN + CaseN`
#> 
#> 
#> 
#> ── 1) Scanning for rows with NA in critical columns ──
#> 
#> 
#> 
#> ✔ No rows contained missing values in CHR, POS, EffectAllele, and OtherAllele
#> 
#> 
#> 
#> ── 2) Scanning for rows with duplications ──
#> 
#> 
#> 
#> ℹ Looking for duplications with columns: CHR, POS, EffectAllele, and OtherAllele
#> 
#> ✔ Found no duplications
#> 
#> 
#> 
#> ── 3) Scanning for indels ──
#> 
#> 
#> 
#> 1. EffectAllele or OtherAllele, character length > 1: A vs AA
#> 
#> 2. EffectAllele or OtherAllele coded as 'D', 'I', or 'R'
#> 
#> ✔ Detected 0 rows as indels
#> 
#> 
#> 
#> ── 4) Validating columns 
#> 
#> ℹ The median value of B is -0.0012, which seems reasonable
#> 
#> 
#> 
#> ── All rows passed validation 
#> 
#> 
#> 
#> ── 5) Adding RSID based on CHR:POS. Adding dbSNP based QC flags 
#> 
#> ℹ Inferring build by matching 10000 rows to GRCh37 and GRCh38
#> 
#> 82 snps matched GRCh38, 9997 for GRCh37, inferring build to be 37
#> ℹ 0 rows had EffectAllele and OtherAllele not matching dbSNP REF/ALT alleles,but which matched after strand-flipping
#> 
#> ! Removed 21 rows with no dbSNP entry or with incompat alleles
#> 
#> /tmp/RtmpRYimIr/file1b6234d4c157/pipeline_info/removed_nodbsnp.parquet
#> 
#> 
#> ── 6) Repairing missings statistics columns if possible ──
#> 
#> 
#> 
#> ℹ Z missing: Calculating Z using the formula:  Z = B / SE
#> 
#> 
#> 
#> ── Finished repair_stats:  
#> 
#> ℹ Added 1 new columns: Z
#> 
#> 
#> 
#> ── Listing final breakdown of removed rows:  
#> 
#> nodbsnp: 21
#> 
#> 
#> 
#> ── Finished tidyGWAS ───────────────────────────────────────────────────────────
#> 
#> ℹ A total of 21 rows were removed
#> 
#> ℹ Total running time: 3.2s
#> 
#> Saving metadata from analysis to /tmp/RtmpRYimIr/file1b6234d4c157/metadata.yaml
#> # A tibble: 99,979 × 21
#>    CHR   POS_37 EffectAllele OtherAllele rowid     B        P     SE  INFO CaseN
#>    <chr>  <int> <chr>        <chr>       <int> <dbl>    <dbl>  <dbl> <dbl> <int>
#>  1 6     2.75e7 C            T           77548 0.212 4.83e-39 0.0162 0.99  53386
#>  2 6     2.81e7 T            C           36396 0.207 4.67e-38 0.016  1     53386
#>  3 6     2.75e7 C            T           77992 0.206 7.92e-38 0.016  0.993 53386
#>  4 6     2.78e7 T            C           67347 0.171 1.45e-32 0.0144 1     53386
#>  5 6     2.63e7 T            C            7577 0.191 8.94e-30 0.0168 0.996 53386
#>  6 6     3.08e7 G            A           26615 0.157 9.77e-28 0.0144 0.954 53386
#>  7 6     2.83e7 C            T           91092 0.116 4.33e-27 0.0108 1.01  53386
#>  8 6     2.64e7 C            T           84820 0.151 9.52e-27 0.0141 1     53386
#>  9 6     3.11e7 C            T           17683 0.143 1.54e-25 0.0137 1     53386
#> 10 6     3.03e7 A            G           71635 0.154 2.58e-25 0.0148 0.992 49683
#> # ℹ 99,969 more rows
#> # ℹ 11 more variables: ControlN <int>, EAF <dbl>, EffectiveN <int>, N <int>,
#> #   POS_38 <int>, RSID <chr>, REF_37 <chr>, REF_38 <chr>, indel <lgl>,
#> #   multi_allelic <lgl>, Z <dbl>
```

#### Output files

tidyGWAS returns all its output in a directory.

``` r
fs::dir_tree(out)
#> /tmp/RtmpRYimIr/file1b6234d4c157
#> ├── metadata.yaml
#> ├── pipeline_info
#> │   ├── removed_nodbsnp.parquet
#> │   └── removed_validate_chr_pos_path.parquet.parquet
#> ├── raw
#> │   └── raw.parquet
#> └── tidyGWAS_hivestyle
#>     ├── CHR=1
#>     │   └── part-0.parquet
#>     ├── CHR=10
#>     │   └── part-0.parquet
#>     ├── CHR=11
#>     │   └── part-0.parquet
#>     ├── CHR=12
#>     │   └── part-0.parquet
#>     ├── CHR=13
#>     │   └── part-0.parquet
#>     ├── CHR=14
#>     │   └── part-0.parquet
#>     ├── CHR=15
#>     │   └── part-0.parquet
#>     ├── CHR=16
#>     │   └── part-0.parquet
#>     ├── CHR=17
#>     │   └── part-0.parquet
#>     ├── CHR=18
#>     │   └── part-0.parquet
#>     ├── CHR=19
#>     │   └── part-0.parquet
#>     ├── CHR=2
#>     │   └── part-0.parquet
#>     ├── CHR=20
#>     │   └── part-0.parquet
#>     ├── CHR=21
#>     │   └── part-0.parquet
#>     ├── CHR=22
#>     │   └── part-0.parquet
#>     ├── CHR=3
#>     │   └── part-0.parquet
#>     ├── CHR=4
#>     │   └── part-0.parquet
#>     ├── CHR=5
#>     │   └── part-0.parquet
#>     ├── CHR=6
#>     │   └── part-0.parquet
#>     ├── CHR=7
#>     │   └── part-0.parquet
#>     ├── CHR=8
#>     │   └── part-0.parquet
#>     └── CHR=9
#>         └── part-0.parquet
```

That’s a lot of files!

1.  Removed rows are stored in `pipeline/removed_*`

2.  `metadata.yaml` contains metadata about the execution

3.  `raw` contains the summary statistics ***before*** any munging was
    done. Useful to reproduce, or to identify why rows were removed.

4.  `tidyGWAS_hivestyle` contains the cleaned summary statistics, in
    something called a [hivestyle
    partition](https://arrow.apache.org/docs/r/articles/dataset.html),
    by default. The motivation for this is detailed further
    [down](#hivestyle-partitioning) in the vignette. If you just want a
    standard csv file, use `output_format="csv"` or
    `output_format="parquet"`

5.  To read in the cleaned summary statistics:

    `df <- arrow::open_dataset(paste0(out, "/tidyGWAS_hivestyle")) |> dplyr::collect()`

#### The tidyGWAS columns

tidyGWAS uses the following column names:

- CHR

- POS

- RSID

- EffectAllele

- OtherAllele

- EAF

- B

- SE

- P

- CaseN

- ControlN

- N

- EffectiveN

- INFO

- Z

**Note**:

If your RSID column is in the format CHR:BP:A1:A2, you can still pass it
as an RSID column. tidyGWAS automatically detects and splits apart the
values.

#### Inputting sample size columns

Often, the sample size column is missing from the summary statistics,
and you provide it manually. tidyGWAS has three arguments that can be
used to manually set the sample size if it’s missing from the original
file:

- `CaseN`

- `ControlN`

- `N`

``` r

cleaned <- tidyGWAS(
  tbl = sumstats, 
  dbsnp_path = dbsnp_path,
  # CaseN, ControlN and N can all be used to set sample size
  CaseN = 400,
  ControlN = 800,
  N = 1200
  )
```

#### Reading in files

If you pass a filepath to tidyGWAS, it will attempt to read in the file
with
[`data.table::fread()`](https://rdatatable.gitlab.io/data.table/reference/fread.html).
If you need to pass custom arguments to fread to make the parsing
correct, you can the arguments through `...`.

``` r
cleaned <- tidyGWAS(
  tbl = "filepath/to/gwas/ondisk.csv",
  # fread has argument `skip`. Let's pretend i want to skip the first 100 rows
  skip = 100.
  dbsnp_path = dbsnp_path,
)
```

## Hivestyle-partitioning

The default output format is a hivestyle format using .parquet files.
See for example
[here](https://arrow.apache.org/docs/r/articles/dataset.html) for some
motivation for this. In essence, this format will significantly speed up
other downstream applications such as meta-analysis, LD querying and
other analyses.

#### Variant identity

TidyGWAS will add the following columns dealing with variant ID, in
addition to saving and validating all valid columns in the input file.

1.  `CHR` Chromosome. The same across both builds

2.  `POS_38`, `POS_37` genomic position on GRCh38 and GRCh37

3.  `RSID` variant ID from dbSNP

4.  `REF_37`, `REF_38` is the reference genome allele on GRCh38 and
    GRCh37

5.  `multi_allelic` is a TRUE/FALSE column that flag rows that were
    multi-allelic IN the summary statistics NOT whether there are
    multiple alleles in dbSNP. (TRUE corresponds to multi allelic).

6.  `rowid` maps each row back to the inputted summary statistics the
    file in `raw`, so that any row can be mapped back to it’s original
    values.

7.  `indel` TRUE/FALSE whether the variant is of type INsertion/DELetion

8.  All other valid columns in the input summary statistics file

9.  If possible
    [`repair_stats()`](https://ararder.github.io/tidyGWAS/reference/repair_stats.md)
    will add statistics columns such as `B, P, SE, Z`, if they are
    missing and possible to repair.

## Parallel computation

tidyGWAS automatically detects the number of cores. In some cases, for
example when running tidyGWAS in a HPC cluster, you might need to
manually set the number of cores, which can be done using the
`OMP_NUM_THREADS` variable. This should not be a larger number of cores
than what you have requested in your HPC job (in the example below, the
“–cpus-per-task” flag)

``` bash
#SBATCH --mem=60gb
#SBATCH --time=24:0:00
#SBATCH --cpus-per-task 8
export OMP_NUM_THREADS=8

outdir=$(pwd)
gwas=$outdir/my_gwas.tsv.gz
dbsnp_files="dbSNP155"
Rscript -e "tidyGWAS(commandArgs(trailingOnly = TRUE)[1],  dbsnp_path = commandArgs(trailingOnly = TRUE)[2],output_dir = commandArgs(trailingOnly = TRUE)[3], logfile=TRUE)" $gwas $dbsnp_files $outdir
```

## Computational cost and memory usage

Memory use and time scales with the size of the summary statistics. From
running tidyGWAS, here’s an estimation. Especially, if the summary
statistics are missing CHR and POS, tidyGWAS will require additional
memory. 1. 10 million rows ~ 5-15gb 2. 40 million rows 10-40gb 3. 60
million rows 65-85gb
