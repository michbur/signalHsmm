% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_hsmm.R
\name{train_hsmm}
\alias{train_hsmm}
\title{Train sighsmm_model object}
\usage{
train_hsmm(train_data, aa_group, max_length = 32,
  region_fun = find_nhc)
}
\arguments{
\item{train_data}{training data.}

\item{aa_group}{method of aggregating amino acids.}

\item{max_length}{maximum length of signal peptide.}

\item{region_fun}{function defining borders of regions (see \code{\link{find_nhc}}).}
}
\value{
object of class \code{sighsmm_model}.
}
\description{
Train sighsmm_model object
}
