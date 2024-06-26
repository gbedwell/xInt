% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/digest.R
\name{digest}
\alias{digest}
\title{Simulate Restriction Digestion}
\usage{
digest(string.list, re.sites, cut.after = 1)
}
\arguments{
\item{string.list}{The output of get_chromosome_seqs() or a comparable list of DNAStrings.}

\item{re.sites}{A character vector of restriction enzyme recognition sequences.}

\item{cut.after}{A single integer that defines the nucleotide after which the DNA is cleaved. Defaults to 1.}
}
\value{
A list of cut site positions across every chromosome for every defined recognition sequence.
}
\description{
Simulate restriction digestion of DNA using arbitrary restriction enzymes. DNA to be digested is input as a list of DNAStrings. Restriction sites are defined in a character vector holding the target recognition sequences.
}
