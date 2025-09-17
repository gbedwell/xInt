#' Get Chromosome Sequences
#'
#' Extracts the nucleotide sequences of each chromosome in a defined genome object
#'
#' @param genome.obj The genome object of interest.
#'
#' @import S4Vectors
#' @import GenomicRanges
#' @import Biostrings
#'
#' @return A list of DNAString objects where each element contains a single chromosome from genome.obj
#'
get_chromosome_seqs <- function(genome.obj){
  seqs <- getSeq(genome.obj, names = names(genome.obj))
  return(seqs)
}
