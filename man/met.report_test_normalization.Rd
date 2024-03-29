% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metabolomics_computation.R
\name{met.report_test_normalization}
\alias{met.report_test_normalization}
\title{Generate a markdown report for the screening of metabolomics data pre-processing methods.}
\usage{
met.report_test_normalization(mSet_list, report.dir = NULL, ...)
}
\arguments{
\item{report.dir}{(Character) Enter the name or location of the folder in which the report files are generated. If \code{NULL} (the default), a new folder "Report_\emph{date-time}" is created in your working directory.}

\item{...}{Further arguments passed from other methods.}

\item{mSetObj}{Enter name of the created mSet object after a full analysis workflow (see \code{\link[VisomX]{met.workflow}}).}
}
\value{
In the specified folder, this function creates a PDF report, an HTML report, and an .RData file including the entire mSet object.
}
\description{
\code{met.report_test_normalization} generates a report of the analysis performed by the \code{met.test_normalization} wrapper function.
}
\author{
Nicolas T. Wirth \email{mail.nicowirth@gmail.com}
Technical University of Denmark
License: GNU GPL (>= 2)
}
