% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow_functions.R
\name{filterAdductGroupsBySize}
\alias{filterAdductGroupsBySize}
\title{Remove features with less than n adducts}
\usage{
filterAdductGroupsBySize(
  x,
  group_size_column = "group_size",
  min_group_size = 3
)
}
\arguments{
\item{x}{MS-DIAL alignment results table}

\item{group_size_column}{name of column (default: "group_size")}

\item{min_group_size}{minimum number of adducts}
}
\value{
an object of the same class as 'x'
}
\description{
Remove features with less than n adducts
}
\examples{
# TBD
}
