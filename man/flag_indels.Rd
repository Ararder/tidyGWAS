% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sumstat_utils.R
\name{flag_indels}
\alias{flag_indels}
\title{Detect "indels" in GWAS summary statistics}
\usage{
flag_indels(tbl)
}
\arguments{
\item{tbl}{a \code{\link[dplyr:reexports]{dplyr::tibble()}} with columns \code{EffectAllele} and \code{OtherAllele}}
}
\value{
a \code{\link[dplyr:reexports]{dplyr::tibble()}} with a TRUE/FALSE column \code{indel} added, where
indel == TRUE corresponds to a row marked as an indel.
}
\description{
Detect "indels" in GWAS summary statistics
}
\examples{
\dontrun{
all_indels <-
  flag_indels(tbl) |>
  dplyr::filter(indels)
}
}
