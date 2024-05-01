#' Make Sequence Logos
#'
#' Generates sequence logos or consensus matrices for integration site datasets.
#'
#'@param site.list The list or GRangesList containing the mapped integration site coordinates.
#'@param seq.len The desired length of the expanded sequences.
#'@param genome.obj The BSgenome object for the genome of interest.
#'@param current.start The position in the target site duplication currently described by the start coordinates in site.list. This is used internally for centering the integration site coordinates.
#'@param tsd The total length of the target site duplication. This is used along with current.start for centering the integration site coordinates.
#'@param return.plot Boolean. Whether or not to return logos for each dataset in site.list instead of the consensus matrix. Defaults to TRUE.
#'
#'@return If return.plot = TRUE, a list of ggplot2 objects containing sequence logos for the mapped integration site datasets in site.list. If return.plot = FALSE, a list of consensus matrices for each integraiton site dataset.
#'
#'@import ggplot2
#'@import ggseqlogo
#'@import GenomicRanges
#'@import Biostrings
#'
#'@export
#'
make_logo <- function( site.list, seq.len = 30, genome.obj, current.start = 3, tsd = 5, return.plot = TRUE ){

  check_sites( site.list )

  if( seq.len %% 2 != 0 ){
    warning("seq.len cannot be odd. Subtracting 1 from seq.len.")
    seq.len <- seq.len - 1
    }

  expanded.sites <- expand_coordinates( site.list = site.list,
                                        seq.len = seq.len,
                                        genome.obj = genome.obj,
                                        current.start = current.start,
                                        tsd = tsd )

  seqs <- getSeq( x = genome.obj,
                  names = expanded.sites,
                  as.character = FALSE )

  mat.list <- lapply( X = seqs,
                      FUN = function(x){
                        consensusMatrix( x, baseOnly=TRUE )
                        }
                      )

  if( isFALSE( return.plot ) ){
    return( mat.list )
  } else{

    center <- ( tsd + 1 ) / 2

    if( tsd %% 2 == 1 ){
      mp <- ( seq.len / 2 )
      zp <- mp - ( center - 1 )
    } else{
      mp <- seq.len / 2
      zp <- mp - ( floor(center) - 1 )
    }

    p.labs <- seq( 1:seq.len ) - zp

    logo.p <- lapply( X = seq_along( mat.list ),
                      FUN = function(x){

                        mat <- mat.list[[x]]

                        ggplot() +
                          geom_logo(data=mat) +
                          theme_bw() +
                          theme(axis.text.y=element_text(size=14),
                                axis.text.x=element_text(size=12),
                                axis.title=element_text(size=16)) +
                          scale_x_continuous(breaks=seq( 1, seq.len, 1 ),
                                             labels=p.labs,
                                             expand=c(0,0)) +
                          scale_y_continuous(expand=c(0,0)) +
                          labs(y="Bits", x="Position", title=names(mat.list)[x])
                        }
                      )

    return( logo.p )
  }
}
