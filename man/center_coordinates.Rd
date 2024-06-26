% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/center_coordinates.R
\name{center_coordinates}
\alias{center_coordinates}
\title{Center Coordinates}
\usage{
center_coordinates(site.list, current.start = 3, tsd = 5, genome.obj)
}
\arguments{
\item{site.list}{The list or GRangesList containing the mapped site coordinates.}

\item{current.start}{The position in the target site duplication currently described by start. This is used for centering the site coordinates.}

\item{tsd}{The total length of the target site duplication. This is used for centering the site coordinates.}

\item{genome.obj}{The genome object of interest.}
}
\value{
A GRanges object containing the centered genomic coordinates.
}
\description{
Centers the input site coordinates relative to the defined target site duplication. In the instance of odd TSDs, returns a 1 bp interval. For even TSDs, a 2 bp interval is returned.
}
