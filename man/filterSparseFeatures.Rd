% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intensity_matrix.R
\name{filterSparseFeatures}
\alias{filterSparseFeatures}
\title{Flag/filter features occurring only in few samples}
\usage{
filterSparseFeatures(
  x,
  sample_group = NULL,
  min_int = 1000,
  min_frac = 0.5,
  flag_column = "passes_sparse_filter",
  filter = FALSE,
  verbose = NULL
)
}
\arguments{
\item{x}{MS-DIAL alignment results table}

\item{sample_group}{numeric or character vector defining sample grouping;
\code{NULL} to use \code{attr(x, "msdial_sam")$class}; "" to apply filter
without sample grouping (equivalent to \code{rep(1, length(samples))}).}

\item{min_int}{minimum intensity that a feature must have, default: 1e3}

\item{min_frac}{minimum fraction of samples (of at least one group if
sample_group is defined accordingly) that contain the feature at an
intensity >= min_int, default: 0.5}

\item{flag_column}{name of flag column, default: "passes_sparse_filter"}

\item{filter}{filter results according to flag, default: FALSE}

\item{verbose}{display messages, defaults to the output of
\code{interactive()}}
}
\value{
an object of the same class as 'x'
}
\description{
Flag/filter features occurring only in few samples
}
