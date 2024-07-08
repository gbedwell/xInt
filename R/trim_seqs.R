#' Trim Fragments
#'
#' Trim the generated random fragments to simulate sequencing reads.
#' Creates paired end fragments of a defined maximum read size and with a defined maximum inner distance.
#'
#'@param fragments The DNAStringSet object containing the random fragment sequences.
#'@param min.width The minimum acceptable read length. Defaults to 14 bp.
#'@param max.distance The maximum desired inner distance. Defaults to 1000 bp.
#'@param max.bp The maximum read length. Defaults to 150 bp.
#'
#'@return A list containing simulated R1 and R2 reads.
#'
#'@import GenomicRanges
#'@import Biostrings
#'
trim_seqs <- function( fragments, min.width = 14, max.distance = 1000, max.bp = 150 ){

  max.width <- max.distance + max.bp*2

  filtered <- fragments[ width( fragments ) >= min.width & width(fragments) <= max.width ]

  left <- subseq( x = filtered, start=1, width = pmin( width( filtered ), max.bp ) )

  right <- subseq(x = filtered,
                  start = pmax( 1, width( filtered ) - max.bp + 1 ), width = pmin( width( filtered ), max.bp ) )

  right <- reverseComplement( right )

  names( left ) <- paste0( "sequence_", seq_along( left ) )

  names( right ) <- paste0( "sequence_", seq_along( right ) )

  return( list( left, right ) )
}
