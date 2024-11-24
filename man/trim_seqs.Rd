% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trim_seqs.R
\name{trim_seqs}
\alias{trim_seqs}
\title{Trim Fragments}
\usage{
trim_seqs(fragments, min.width = 10, max.width = 1000, max.bp = 150)
}
\arguments{
\item{fragments}{The DNAStringSet object containing the random fragment sequences.}

\item{min.width}{The minimum acceptable fragment length (insert size). Defaults to 10.}

\item{max.width}{The maximum acceptable fragment length (insert size). Defaults to 1000.}

\item{max.bp}{The maximum read length. Defaults to 150 bp.}
}
\value{
A list containing simulated R1 and R2 reads.
}
\description{
Trim the generated random fragments to simulate sequencing reads.
Creates paired end fragments of a defined maximum read size and with a defined maximum inner distance.
}
