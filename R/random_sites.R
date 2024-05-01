#' Generate Random Sites
#'
#' Generates randomly positioned integration sites within the genome of interest.
#'
#'@param n.sites The number of sites to generate.
#'@param genome.obj The BSgenome object of interest.
#'
#'@return A GRanges object containing the integration site positions.
#'
#'@import GenomeInfoDb
#'@import GenomicRanges
#'@importFrom IRanges IRanges
#'
random_sites <- function( n.sites, genome.obj ){

  chr.names <- seqlevels( genome.obj )
  chr.lengths <- seqlengths( genome.obj )
  chr.weights <- chr.lengths / sum( chr.lengths )
  chr.vec <- sample( x = chr.names, size = n.sites, prob = chr.weights, replace = TRUE )
  site.pos <- vapply( X = chr.vec,
                      FUN = function(x) {
                        cl <- chr.lengths[[x]]
                        sample( x = seq_len( cl ), size = 1 )
                        },
                      FUN.VALUE = numeric(1)
                      )
  site.strand <- sample( x = c( "+", "-" ), size = n.sites, replace = TRUE )
  sitesG <- GRanges( seqnames = chr.vec,
                     ranges = IRanges( start = site.pos, width = 1 ),
                     strand = site.strand )
  return( sitesG )

}
