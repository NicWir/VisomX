% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteomics_computation.R
\name{manual_impute}
\alias{manual_impute}
\title{Imputation by random draws from a manually defined distribution}
\usage{
manual_impute(se, scale = 0.3, shift = 1.8)
}
\arguments{
\item{se}{SummarizedExperiment, Proteomics data (output from make_se() or make_se_parse()). It is adviced to first remove proteins with too many missing values using filter_missval() and normalize the data using normalize_vsn().}

\item{scale}{Numeric(1), Sets the width of the distribution relative to the standard deviation of the original distribution.}

\item{shift}{Numeric(1), Sets the left-shift of the distribution (in standard deviations) from the median of the original distribution.}
}
\value{
A SummarizedExperiment object with imputed values
}
\description{
This function imputes missing values in a proteomics dataset by random draws from a manually defined distribution.
}
