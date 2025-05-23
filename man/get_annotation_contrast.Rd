% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_misc_utils.R
\name{get_annotation_contrast}
\alias{get_annotation_contrast}
\title{Subset and annotate a DESeqDataSet for contrast-specific heatmaps}
\usage{
get_annotation_contrast(dds, indicate, contrast = contrast_samples)
}
\arguments{
\item{dds}{A \code{DESeqDataSet} or \code{SummarizedExperiment}.}

\item{indicate}{Character vector of column names in \code{colData(dds)} to use for annotation.}

\item{contrast}{Character vector of condition patterns to retain in the annotation.}
}
\value{
A \code{ComplexHeatmap::HeatmapAnnotation} object.
}
\description{
Subset and annotate a DESeqDataSet for contrast-specific heatmaps
}
