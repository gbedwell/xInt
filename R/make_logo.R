#' Make Sequence Logos
#'
#' Generates sequence logos for each dataset in site.list.
#'
#'@param site.list The list or GRangesList containing the mapped site coordinates.
#'@param seq.len The length of the expanded sequences.
#'@param genome.obj The genome object of interest.
#'@param current.start The position in the target site duplication currently described by start. This is used for centering the site coordinates.
#'@param tsd The total length of the target site duplication. This is used for centering the site coordinates.
#'@param return.plot Boolean. Whether or not to return diagnostic plots for each dataset in site.list instead of the consensus matrix. Defaults to TRUE.
#'
#'@export
#'
make_logo <- function( site.list, seq.len = 50, genome.obj, current.start = 3, tsd = 5, return.plot = TRUE ){

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
                  names = as( expanded.sites, "GRangesList" ),
                  as.character = FALSE )

  mat.list <- lapply( X = seqs,
                      FUN = function(x){
                        Biostrings::consensusMatrix( x, baseOnly=TRUE)
                        }
                      )

  if( isFALSE( return.plot ) ){
    return( mat.list )
  } else{

    center <- ( tsd + 1 ) / 2

    if( tsd %% 2 == 1 ){
      mp <- ( seq.len / 2 ) + 1
      zp <- mp - ( center )
    } else{
      mp <- seq.len / 2
      zp <- mp - ( floor(center) )
    }

    p.labs <- seq( 1:seq.len ) - zp

    logo.p <- lapply( X = mat.list,
                      FUN = function(x){

                        require(ggplot2)
                        require(ggseqlogo)

                        ggplot() +
                          geom_logo(data=x) +
                          theme_bw() +
                          theme(axis.text.y=element_text(size=14),
                                axis.text.x=element_text(size=12),
                                axis.title=element_text(size=16)) +
                          scale_x_continuous(breaks=seq( 1, seq.len, 1 ),
                                             labels=p.labs,
                                             expand=c(0,0)) +
                          scale_y_continuous(expand=c(0,0)) +
                          labs(y="Bits", x="Position")
                        }
                      )

    return( logo.p )
  }
}
