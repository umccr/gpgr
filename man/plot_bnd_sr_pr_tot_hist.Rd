% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sv.R
\name{plot_bnd_sr_pr_tot_hist}
\alias{plot_bnd_sr_pr_tot_hist}
\title{Histogram for SR, PR and SR + PR for BNDs}
\usage{
plot_bnd_sr_pr_tot_hist(
  d,
  title = "SR, PR and SR + PR histogram for BNDs",
  subtitle = "Values of 0 (NA) are not shown."
)
}
\arguments{
\item{d}{A data.frame with an SR_PR_alt column.}

\item{title}{Main title of plot.}

\item{subtitle}{Subtitle of plot.}
}
\value{
A ggplot2 plot object.
}
\description{
Plots histograms for the number of split reads (\code{SR}), paired end reads (\code{PR}), and their
sum (\code{tot}) across all BNDs. Observations where the SR or PR value is 0 (NA) are not shown.
}
\examples{
x <- system.file("extdata/sash/sv.prioritised.tsv", package = "gpgr")
d <- process_sv(x)$map
plot_bnd_sr_pr_tot_hist(d, "a title")
}
