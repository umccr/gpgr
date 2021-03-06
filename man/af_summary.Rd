% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/umccrise.R
\name{af_summary}
\alias{af_summary}
\title{umccrise AF summary}
\usage{
af_summary(af_global_file, af_keygenes_file)
}
\arguments{
\item{af_global_file}{Path to 'global' AF file.}

\item{af_keygenes_file}{Path to 'keygenes' AF file.}
}
\value{
A list containing a gt table and a plot summarising AFs.
}
\description{
Get summary table and plot for SNV allele frequencies output by umccrise.
}
