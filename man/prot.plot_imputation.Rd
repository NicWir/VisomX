% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteomics_plots.R
\name{prot.plot_imputation}
\alias{prot.plot_imputation}
\title{Plot imputation results}
\usage{
prot.plot_imputation(se, ..., plot = TRUE, basesize = 12, export = FALSE)
}
\arguments{
\item{se}{SummarizedExperiment object}

\item{...}{Additional SummarizedExperiment objects. Can be input as named argument, for example: \code{"se2" = se2}}

\item{plot}{logical, whether to plot the results or not (default: TRUE)}

\item{basesize}{numeric, base font size (default: 12)}

\item{export}{logical, whether to export the results as PNG and PDF or not (default: TRUE)}
}
\value{
(invisible) a ggplot object
}
\description{
Plot the results of imputation as a density plot. It is recommended to log-transform the data before using this function.
}
