% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metabolomics_plots.R
\name{met.plot_PLS3DScore}
\alias{met.plot_PLS3DScore}
\title{Create 3D PLS-DA Score plot}
\usage{
met.plot_PLS3DScore(
  mSetObj = NA,
  imgName = "PLS3DScore",
  format = "json",
  inx1 = 1,
  inx2 = 2,
  inx3 = 3,
  export = F
)
}
\arguments{
\item{mSetObj}{Input name of the created mSet object.
Data container after partial least squares-discriminant analysis (\code{\link[VisomX]{met.PLSR.Anal}} and \code{\link[VisomX]{met.PLSDA.CV}}).}

\item{imgName}{(Character) Enter a name for the image file (if \code{export = TRUE}).}

\item{format}{format.}

\item{inx1}{(Numeric) Indicate the number of the principal component for the x-axis of the score plot.}

\item{inx2}{(Numeric) Indicate the number of the principal component for the y-axis of the score plot.}

\item{inx3}{(Numeric) Indicate the number of the principal component for the z-axis of the score plot.}

\item{export}{(Logical, \code{TRUE} or \code{FALSE}) Shall the plot be exported as file with the chosen format?}
}
\value{
The input mSet object with added 3D scatter plot. The plot can be retrieved from within R by executing \code{met.print_PLS3DScore(mSetObj$imgSet$pls.score3d.plot)}.
}
\description{
\code{met.plot_PLS3DScore} visualizes clusters of samples based on their similarity in principal component analysis.
}
\author{
Nicolas T. Wirth \email{mail.nicowirth@gmail.com}
}
