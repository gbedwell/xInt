#' Expand Coordinates
#'
#' Expands site coordinates by a defined amount.
#'
#'@param site.list The list or GRangesList containing the mapped site coordinates.
#'@param seq.len The length of the returned sequence.
#'@param genome.obj The genome object of interest.
#'@param current.start The position in the target site duplication currently described by start.
#'This is used for centering the site coordinates.
#'@param tsd The total length of the target site duplication. This is used for centering the site coordinates.
#'
#'@return A GRanges object containing the expanded genomic coordinates.
#'
#'@import methods
#'@import GenomicRanges
#'
expand_coordinates <- function( site.list, seq.len = 50, genome.obj, current.start = 1, tsd = 5 ){

  if( seq.len %% 2 != 0 ){
    warning("seq.len cannot be odd. Subtracting 1 from seq.len.")
    seq.len <- seq.len - 1
  }

  cc <- center_coordinates( site.list = site.list, current = current.start, tsd = tsd, genome.obj = genome.obj )

  ex <- lapply( X = cc,
                FUN = function(x){
                  # Anchors the expanded sequence on the central coordinate
                  # The effect of this is to shorten seq.len by 1 when the central base is integer-valued (tsd is odd).
                  # This is different from GenomicRanges::resize(), which will add 1 more base to left side
                  # of odd-width coordinates.
                  start(x) <- start(x) - ((seq.len / 2) - 1)
                  end(x) <- end(x) + ((seq.len / 2) - 1)
                  outs <- bound_check( fragments = x,
                                       genome.obj = genome.obj,
                                       include.lower = TRUE )
                  if( length(outs) != 0 ){
                    return(x[-outs])
                  } else{
                    return(x)
                  }
                }
  )

  if( is( site.list, "GRangesList" ) ){
    ex <- as( ex, "GRangesList" )
  }

  return( ex )
}
