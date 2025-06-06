% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectrum_functions.r
\name{cleanSpec}
\alias{cleanSpec}
\title{Clean mass spectrum}
\usage{
cleanSpec(x, relthr = 0, absthr = 0, dmz = NULL, minmz = NULL, maxmz = NULL)
}
\arguments{
\item{x}{matrix or data.frame}

\item{relthr}{relative intensity threshold}

\item{absthr}{absolute intensity threshold}

\item{dmz}{delta m/z}

\item{minmz}{minimum m/z}

\item{maxmz}{maximum m/z}
}
\value{
matrix or data.frame
}
\description{
Clean mass spectrum
}
