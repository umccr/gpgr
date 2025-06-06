% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sv.R
\name{plot_bnd_sf_df_tot_lines}
\alias{plot_bnd_sf_df_tot_lines}
\title{Line plot for SF, DF and SF + DF for BNDs}
\usage{
plot_bnd_sf_df_tot_lines(
  d,
  title = "SF, DF and SF + DF  line plot for BNDs",
  subtitle = "Events are sorted by decreasing total values."
)
}
\arguments{
\item{d}{A data.frame with SF_alt, DF_alt, and SF_DF_alt columns.}

\item{title}{Main title of plot.}

\item{subtitle}{Subtitle of plot.}
}
\value{
A ggplot2 plot object.
}
\description{
Plots the split fragments (\code{SF}), discordant fragments (\code{DF}), and
the total fragments (SF+DF) across all BNDs, sorted by total value.
}
\examples{
x <- system.file("extdata/sash/sv.prioritised.tsv", package = "gpgr")
d <- process_sv(x)$map
plot_bnd_sf_df_tot_lines(d)
}
