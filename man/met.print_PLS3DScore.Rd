% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metabolomics_plots.R
\name{met.print_PLS3DScore}
\alias{met.print_PLS3DScore}
\title{Plot a generated 3D PLS-DA scores plot in RStudio}
\usage{
met.print_PLS3DScore(mSetObj = NA)
}
\arguments{
\item{mSetObj}{Input name of the created mSet object.
Data container after first after \code{\link[VisomX]{met.plot_PLS3DScore}}.}
}
\value{
A 3D scatter plot. To close the graphics device, execute \code{rgl::close3d()}.
}
\description{
\code{met.print_PLS3DScore} opens a new graphics device to show the generated 3D PLS-DA score plot.
}
\author{
Nicolas T. Wirth \email{mail.nicowirth@gmail.com}
}
