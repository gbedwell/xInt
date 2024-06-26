% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_metadata.R
\name{add_metadata}
\alias{add_metadata}
\title{Add Metadata}
\usage{
add_metadata(
  xint.obj = NULL,
  virus,
  cell.type,
  poi,
  num.cells,
  dpi,
  moi,
  frag.method,
  seq.platform,
  seq.type,
  read.length,
  date
)
}
\arguments{
\item{xint.obj}{The xInt object of interest.}

\item{virus}{The type of virus used in the infection.}

\item{cell.type}{The type of cell infected.}

\item{poi}{The protein-of-interest in transfected cells.}

\item{num.cells}{The number of cells infected.}

\item{dpi}{The days post infection.}

\item{moi}{The multiplicity of infection.}

\item{frag.method}{The method of genome fragmentation.}

\item{seq.platform}{The sequencing platform used.}

\item{seq.type}{The type of sequencing (e.g., paired- or single-end).}

\item{read.length}{The maximum read length.}

\item{date}{The date of sequencing.}
}
\value{
The input xInt object with metadata attached.
}
\description{
Adds basic experiment metadata to an xInt object. Intended to be used to attach important information (e.g., virus type, cell type, genome fragmentation method, etc.) directly to the related xInt object. Any/all metadata fields can be left blank, if desired.
}
