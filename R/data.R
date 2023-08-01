#' PGC3 schizophrenia GWAS
#'
#' A subset of data from the Psychiatics Genonomics Consortium GWAS on schizophrenia (european subset)
#'
#' @format ## `test_file`
#' A data frame with 100,000 rows and 12 columns:
#' \describe{
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic position}
#'   \item{RSID}{rsID from dbSNP}
#'   \item{EffecttAllele}{The allele that correspods to the effect, B}
#'   \item{OtherAllelle}{Non-effect allele}
#'   \item{B}{Effect, Beta}
#'   \item{SE}{Standard error of B}
#'   \item{EAF}{EffectAllele Frequency, frequency of EffectAllele}
#'   \item{INFO}{Imputation accuracy}
#'   \item{P}{Pvalue}
#'   \item{CaseN}{Number of cases}
#'   \item{ControlN}{Number of controls}
#'   ...
#' }
#' @source <https://figshare.com/ndownloader/files/34517828>
"test_file"


#' RsMergeArch
#'
#' A small subset of RsMergeArch dbsNP 151, with updated format. Consists of the 175 rows
#' where test_file has an RSID that has been merged
#'
#' @format ## `rs_merge_arch`
#' A data frame with 172 rows and 2 columns:
#' \describe{
#'   \item{RSID}{old, retracted RSID}
#'   \item{new_RSID}{new RSID, the one which remains after merge with RSID}
#' }
#' @source <https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/database/organism_data/RsMergeArch.bcp.gz>
"rs_merge_arch"
