% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_misc_utils.R
\name{.onAttach}
\alias{.onAttach}
\title{Run upon attaching package VisomX}
\usage{
.onAttach(libname, pkgname)
}
\arguments{
\item{libname}{library name}

\item{pkgname}{package name}
}
\description{
Changes debug option for package \code{rgl} to avoid Rstudio crashing upon attaching it and prints welcome message
}
