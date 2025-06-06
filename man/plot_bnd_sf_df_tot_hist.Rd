% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sv.R
\name{plot_bnd_sf_df_tot_hist}
\alias{plot_bnd_sf_df_tot_hist}
\title{Histogram for SF, DF and SF_DF_sum for BNDs}
\usage{
plot_bnd_sf_df_tot_hist(
  d,
  title = "SF, DF and SF + DF histogram for BNDs",
  subtitle = "Values of 0 (NA) are not shown."
)
}
\arguments{
\item{d}{A data.frame with SF_alt, DF_alt, and SF_DF_sum columns.}

\item{title}{Main title of plot.}

\item{subtitle}{Subtitle of plot.}
}
\value{
A ggplot2 plot object.
}
\description{
Plots histograms for split fragments (\code{SF}), discordant fragments (\code{DF}), and
the total fragments (SF+DF) across all BNDs. Observations where values are 0 (NA) are not shown.
}
\examples{
x <- system.file("extdata/sash/sv.prioritised.tsv", package = "gpgr")
d <- process_sv(x)$map
plot_bnd_sf_df_tot_hist(d, "a title")
}
