#' Pool Replicate Datasets
#'
#' Given a list of integration site datasets containing replicates from various conditions,
#'pool the replicates into a single, larger dataset while keeping different conditions separate.
#'
#'@param site.list A list of GRanges objects or a GRangesList containing integration site coordinates from various samples.
#'@param group.index.list A list of index vectors enumerating which site.list elements correspond to which conditions.
#'Each group.index.list element should correspond to a unique condition.
#'
#'@return A list of GRanges objects or a GRangesList of the combined coordinates.
#'
#'@examples
#'data(sites)
#'collapse_sites(site.list = sites,
#'               group.index.list = list(A = 1:4,
#'                                       B = 5:9,
#'                                       C = 10)
#'
#'@import methods
#'@import GenomicRanges
#'
#'@export
#'
collapse_sites <- function( site.list, group.index.list ){

  if( !is.list( group.index.list ) ){
    stop( "Group index values must be given as a list.",
          call. = FALSE )
  }

  if( !validObject( site.list ) ){
    stop( "site.list is not a valid SiteListObject.",
          call. = FALSE )
  }

  ll <- lapply( X = group.index.list,
                FUN = function(x){
                  ind <- x
                  out <- unlist( as( site.list[ ind ], "GRangesList" ) )
                }
  )

  if( is.null( names( ll ) ) ){
    warning( "No group names were provided in
             'group.index.list. The returned list
             will therefore not be named.",
             call. = FALSE )
  }

  if( is( site.list, "GRangesList" ) ){
    ll <- as( ll, "GRangesList" )
  }

  return( ll )
}
