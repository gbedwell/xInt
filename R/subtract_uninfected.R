#' Subtract Uninfected Sites
#'
#' Subtract integration sites spuriously obtained from uninfected control samples.
#' Uninfected sites are collapsed into a single
#'
#'@param site.list A list of GRanges objects or a GRangesList containing integration site datasets.
#'@param uninfected.datasets A character vector containing the names of uninfected datasets in site.list
#'or a numeric vector containing the indeces of uninfected datasets in site.list.
#'These datasets are removed in the output when return.uninfected = FALSE.
#'@param return.uninfected Boolean. Whether or not to return the combined uninfected site coordinates.
#'
#'@return A list of GRanges objects or a GRangesList containing the subtracted integration site datasets.
#'
#'@import methods
#'@import GenomicRanges
#'
#'@export
#'
subtract_uninfected <- function( site.list,
                                 uninfected.datasets,
                                 return.uninfected = FALSE ){

  if( !validObject( site.list ) ){
    stop( "site.list is not a valid SiteListObject.",
          call. = FALSE )
  }

  if( !is.character( uninfected.datasets ) ){
    if( is.numeric( uninfected.datasets ) ){
      uninfected.datasets <- names( site.list )[ uninfected.datasets ]
    } else{
      stop( "uninfected.datasets must be a numeric or character vector.",
            call. = FALSE )
    }
  }

  spur <- site.list[ which( names(site.list) %in% uninfected.datasets ) ]

  if( is.list( site.list ) ){
    spur <- unlist( as( spur, "GRangesList" ) )
  } else{
    spur <- unlist( site.list )
  }


  if( isTRUE( return.uninfected ) ){
    return( spur )
  } else{
    sub.sites <- lapply( X = site.list[ -which( names(site.list) %in% uninfected.datasets ) ],
                         FUN = function(x){
                           ss <- subtract( x = x, y = spur, minoverlap = 1L, ignore.strand = FALSE )
                           unlist( as( ss, "GRangesList" ) )
                           }
                         )

    if( is( site.list, "GRangesList" ) ){
      sub.sites <- as( sub.sites, "GRangesList" )
    }

    return(sub.sites)
  }
}
