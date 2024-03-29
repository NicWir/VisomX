% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transcriptomics_computation.R
\name{rna.filter_missing}
\alias{rna.filter_missing}
\title{Filter genes based on missing values}
\usage{
rna.filter_missing(
  se,
  type = c("complete", "condition", "fraction", NULL),
  thr = NULL,
  min = NULL
)
}
\arguments{
\item{se}{\code{SummarizedExperiment} object, transcriptomics data parsed with \code{\link{rna.read_data}}.}

\item{type}{(Character string) "complete", "condition" or "fraction", Sets the type of filtering applied. "complete" will only keep genes with valid values in all samples. "condition" will keep genes that have a maximum of \code{thr} missing values in at least one condition. "fraction" will keep genes that have a \code{min} fraction of valid values in all samples.}

\item{thr}{(Integer) Sets the threshold for the allowed number of missing values in at least one condition if \code{type = "condition"}. In other words: "keep genes that have a maximum of 'thr' missing values in at least one condition."}

\item{min}{(Numeric) Sets the threshold for the minimum fraction of valid values allowed for any gene if \code{type = "fraction"}.}
}
\value{
A filtered SummarizedExperiment object.
}
\description{
\code{rna.filter_missing} filters a transcriptomics dataset based on missing values. Different types of filtering can be applied, which range from only keeping genes without missing values to keeping genes with a certain percent valid values in all samples or keeping genes that are complete in at least one condition.
}
