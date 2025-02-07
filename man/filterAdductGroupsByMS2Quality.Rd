% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow_functions.R
\name{filterAdductGroupsByMS2Quality}
\alias{filterAdductGroupsByMS2Quality}
\title{Remove features with bad MS/MS spectra}
\usage{
filterAdductGroupsByMS2Quality(
  x,
  feat_group_column = "feat_group",
  ms2_column = "ms2",
  min_tic = 1000,
  min_ions = 3,
  keep_all_adducts = FALSE
)
}
\arguments{
\item{x}{MS-DIAL alignment results table}

\item{feat_group_column}{name of column (default: "feat_group")}

\item{ms2_column}{name of column (default: "ms2")}

\item{min_tic}{minimum summed intensity of MS/MS spectrum}

\item{min_ions}{minimum number of ions of MS/MS spectrum}

\item{keep_all_adducts}{keep all or only adduct with best MS/MS spectrum
within each feature group}
}
\value{
an object of the same class as 'x'
}
\description{
Remove features with bad MS/MS spectra
}
\examples{
# TBD
}
