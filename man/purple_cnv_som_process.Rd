% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/purple.R
\name{purple_cnv_som_process}
\alias{purple_cnv_som_process}
\title{Process PURPLE CNV Somatic File for UMCCRISE}
\usage{
purple_cnv_som_process(x)
}
\arguments{
\item{x}{Path to \code{purple.cnv.somatic.tsv} file.}
}
\value{
List with two elements:
\itemize{
\item \code{tab}: Tibble with more condensed columns.
\item \code{descr}: Description of tibble columns.
}
}
\description{
Processes the \code{purple.cnv.somatic.tsv} file and selects columns of interest.
}
\examples{
x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
(pp <- purple_cnv_som_process(x))
}
