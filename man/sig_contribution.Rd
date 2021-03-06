% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutationalpatterns.R
\name{sig_contribution}
\alias{sig_contribution}
\title{Calculate COSMIC Signature Contribution}
\usage{
sig_contribution(mut_mat, signatures)
}
\arguments{
\item{mut_mat}{Mutation count matrix (dimensions: m rows (mutation types) X 1 column (sample)).}

\item{signatures}{Signature matrix (dimensions: m rows (mutation types) X n columns (signatures))}
}
\value{
A list with the COSMIC 2015 and 2020 signature contributions to the
sample's signature.
}
\description{
Finds the linear combination of COSMIC (2015 and 2020) mutation signatures that
most closely reconstructs the SNV mutation matrix by solving the
nonnegative least-squares constraints problem.
}
