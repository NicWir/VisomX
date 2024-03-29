% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteomics_computation.R
\name{prot.normalize_vsn}
\alias{prot.normalize_vsn}
\title{Normalize the data via variance stabilization normalization}
\usage{
prot.normalize_vsn(se, plot = TRUE, export = TRUE)
}
\arguments{
\item{se}{A SummarizedExperiment object}

\item{plot}{Logical. If TRUE, plots the meanSdPlot}

\item{export}{Logical. If TRUE, exports the meanSdPlot in pdf and png}
}
\value{
A SummarizedExperiment object
}
\description{
Normalize the data via variance stabilization normalization
}
