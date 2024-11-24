% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expand_coordinates.R
\name{expand_coordinates}
\alias{expand_coordinates}
\title{Expand Coordinates}
\usage{
expand_coordinates(
  site.list,
  seq.len = 50,
  genome.obj,
  current.start = 1,
  tsd = 5
)
}
\arguments{
\item{site.list}{The list or GRangesList containing the mapped site coordinates.}

\item{seq.len}{The length of the returned sequence.}

\item{genome.obj}{The genome object of interest. This is used for checking the validity of
the centered and expanded coordinates against the genome boundaries.
Coordinates that violate genome boundaries are removed.
This is strictly enforced with expanded coordinates because sequence logos, etc.
require sequences of the same length.}

\item{current.start}{The position in the target site duplication currently described by start.
This is used for centering the site coordinates.}

\item{tsd}{The total length of the target site duplication. This is used for centering the site coordinates.}
}
\value{
A GRanges object containing the expanded genomic coordinates.
}
\description{
Expands site coordinates by a defined amount.
}
