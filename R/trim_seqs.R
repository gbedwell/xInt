#' Trim Fragments
#'
#' Trims the fragments to simulate the sequencing run (e.g. 150 bp reads). Creates paired end fragments of a defined maximum size and with a defined inner distance.
#'
#'@param fragments The fragments generated with BSgenome::getSeq().
#'@param min.width The minimum desired fragment/read length. Default 14 bp.
#'@param max.distance The maximum desired inner distance. Default 1000 bp.
#'@param max.bp The maximum read length. Default 150 bp.
#'
#'@import Biostrings
#'
#'@export
#'
trim_seqs <- function( fragments, min.width = 40, max.distance = 1000, max.bp = 150 ){
  max.width <- max.distance + max.bp*2
  filtered <- fragments[ width( fragments ) >= min.width & width(fragments) <= max.width ]
  left <- subseq( x = filtered, start=1, width = pmin( width( filtered ), 150 ) )
  right <- subseq(x = filtered,
                  start = pmax( 1, width( filtered ) - 150 + 1 ), width = pmin( width( filtered ), 150 ) )
  right <- reverseComplement( right )
  names( left ) <- paste0( "sequence_", seq_along( left ) )
  names( right ) <- paste0( "sequence_", seq_along( right ) )
  return( list( left, right ) )
}
