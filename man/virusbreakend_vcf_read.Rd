% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oncoviral.R
\name{virusbreakend_vcf_read}
\alias{virusbreakend_vcf_read}
\title{Read VIRUSBreakend VCF File}
\usage{
virusbreakend_vcf_read(x)
}
\arguments{
\item{x}{Path to \code{virusbreakend.vcf} file.}
}
\value{
List with two elements:
\itemize{
\item \code{tab}: Tibble containing selected data.
\item \code{descr}: Description of tibble columns.
}
}
\description{
Reads the \code{virusbreakend.vcf} file and selects data to present.
}
\examples{
x <- system.file("extdata/virusbreakend/virusbreakend.vcf", package = "gpgr")
(vb <- virusbreakend_vcf_read(x))
}
