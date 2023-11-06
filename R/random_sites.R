#' Generate Random Sites
#'
#' Generates random integration sites within the defined genome
#'
#'@param n.sites The number of sites to generate.
#'@param genome.obj The name of the genome object of interest.
#'
#'@import GenomeInfoDb
#'@import GenomicRanges
#'
#'@export
#'
random_sites <- function( n.sites, genome.obj ){

  chr.names <- seqlevels( genome.obj )
  chr.lengths <- seqlengths( genome.obj )
  chr.weights <- chr.lengths / sum( chr.lengths )
  chr.vec <- sample( chr.names, n.sites, prob = chr.weights, replace = TRUE )
  site.pos <- sapply( X = chr.vec,
                      FUN = function(x) {
                        sample(x = 1:chr.lengths[[x]], size = 1 )
                        }
                      )
  site.strand <- sample( x = c( "+", "-" ), size = n.sites, replace = TRUE )
  sitesG <- GRanges( seqnames = chr.vec,
                     ranges = IRanges( start=site.pos, width=1 ),
                     strand = site.strand )
  return( sitesG )

}
