#' Genome Boundary Check
#'
#' Check that ranges do not violate genome/chromosome boundaries.
#'
#' @param fragments The GRanges object being tested.
#' @param genome.obj The BSgenome object of interest.
#' @param include.lower Boolean. Whether or not the check lower bound validity.
#'
#' @return A vector of index values for terms in violation of genomic boundaries.
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#'
bound_check <- function(fragments, genome.obj, include.lower = FALSE){
  frag.seqnames <- as.character(GenomicRanges::seqnames(fragments))
  genome.seqlengths <- GenomeInfoDb::seqlengths(genome.obj)
  if(isTRUE(include.lower)){
    which(GenomicRanges::end(fragments) > genome.seqlengths[c(frag.seqnames)] | start(fragments) < 1)
  } else{
    which(GenomicRanges::end(fragments) > genome.seqlengths[c(frag.seqnames)])
  }
}
