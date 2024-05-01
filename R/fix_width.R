#' Fix Coordinate Widths
#'
#' Shift integration site positions to the desired position in the target site duplication and enforce that their width is 1.
#'
#'@param site.list The list of GRanges objects or GRangesList containing the mapped integration site coordinates.
#'@param current.pos The position in the target site duplication described by the current start position. Defaults to 1.
#'@param target.pos The desired position of the shifted start position in the target site duplication. Defaults to 1.
#'@param genome.obj The BSgenome object of interest.
#'
#'@return A list of GRanges objects or a GRangesList containing coordinates with width = 1.
#'
#'@import methods
#'@import GenomicRanges
#'
#'@export
#'
fix_width <- function( site.list, current.pos = 1, target.pos = 1, genome.obj ){

  delta <- target.pos - current.pos

  if( delta %% 1 != 0 ){
    stop( "current.pos and target.pos must both be integers.",
          call. = FALSE )
  }

  sites <- lapply(
    X = site.list,
    FUN = function(x){
      plus <- x[ strand(x) == "+" ]
      minus <- x[ strand(x) == "-" ]
      start( plus ) <- start( plus ) + delta
      end( plus ) <- start( plus )
      end( minus ) <- end( minus ) - delta
      start( minus ) <- end( minus )

      tmp <- sort( c( plus, minus ), ignore.strand = TRUE )

      outs <- bound_check( fragments = tmp,
                           genome.obj = genome.obj,
                           include.lower = TRUE )

      if( length(outs) != 0 ){

        warning( "Shifted coordinates " , paste(outs, collapse=", "),
                 " are out of bounds. These coordinates will be removed.",
                 call. = FALSE )

        tmp <- tmp[-outs]
      }
      return( tmp )
    }
  )

  if( is( site.list, "GRangesList" ) ){
    sites <- as( sites, "GRangesList" )
  }

  return( sites )
}
