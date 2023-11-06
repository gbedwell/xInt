#' Genome Boundary Check
#'
#' Check that simulated random fragment ranges do not violate genome/chromosome boundaries. Not intended to be run alone. Called within make_fragments().
#'
#'@param fragments The random fragments being tested.
#'@param genome.obj The name of the genome object of interest.
#'
#'@returns A vector of index values for terms in violation of genomic boundaries.
#'
#'@importFrom GenomicRanges "seqnames"
#'@importFrom GenomicRanges "end"
#'@importFrom GenomeInfoDb "seqlengths"
#'
#'@examples
#'\dontrun{
#'bound_check(fragments = frag.object, genome.obj = genome.object)
#'}
#'
#'@export
#'
bound_check <- function( fragments, genome.obj ){
  frag.seqnames <- as.character( GenomicRanges::seqnames( fragments ) )
  genome.seqlengths <- GenomeInfoDb::seqlengths( genome.obj )
  which( GenomicRanges::end( fragments ) > genome.seqlengths[ c( frag.seqnames ) ] )
}
