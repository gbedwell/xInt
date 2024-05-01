#' Get Chromosome Sequences
#'
#' Extracts the nucleotide sequences of each chromosome in a defined genome object
#'
#'@param genome.obj The genome object of interest.
#'
#'@import S4Vectors
#'@import GenomicRanges
#'@import Biostrings
#'
#'@return A list of DNAString objects where each element contains a single chromosome from genome.obj
#'
get_chromosome_seqs <- function( genome.obj ){
  seqs <- lapply( X = seqlevels( genome.obj ),
                  FUN = function(x){
                    string <- DNAString( genome.obj[[ x ]] )
                    metadata( string )$name <- x
                    return( string )
                    } )
  names( seqs ) <- seqnames( genome.obj )
  return( seqs )
}
