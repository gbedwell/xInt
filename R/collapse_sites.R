#' Pool Replicate Datasets
#'
#' Given a list of integration site datasets containing replicates from various conditions, pool the replicates into a single, larger dataset while keeping different conditions separate.
#'
#'@param site.list A list of GRanges objects or a GRangesList containing integration site coordinates from various samples.
#'@param group.index.list A list of index vectors enumerating which site.list elements correspond to which conditions. Each group.index.list element should correspond to a unique condition.
#'@param group.names A character vector of names corresponding to the represented conditions. The order should match the order given in group.index.list.
#'
#'@return A list of GRanges objects or a GRangesList of the combined coordinates.
#'
#'@import methods
#'@import GenomicRanges
#'
#'@export
#'
collapse_sites <- function( site.list, group.index.list, group.names = NULL ){

  if( !is.list( group.index.list ) ){
    stop( "Group index values must be given as a list.",
          call. = FALSE )
  }

  check_sites( site.list )

  ll <- lapply( X = seq_along( group.names ),
                FUN = function(x){
                  ind <- group.index.list[[x]]
                  out <- unlist( as( site.list[ ind ], "GRangesList" ) )
                }
  )

  if( !is.null( group.names ) ){
    if( length( group.index.list ) != length( group.names ) ){
      stop( "The number of groups should equal the given group names.",
            call. = FALSE )
    }
    names( ll ) <- group.names
  }

  if( is.list( site.list ) ){
    return( ll )
  } else{
    ll <- as( ll, "GRangesList" )
    return(ll)
  }
}
