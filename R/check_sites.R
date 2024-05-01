#' Check Sites for Compatibility
#'
#' Verify that imported IS datasets meet the requirements for downstream xInt analyses.
#'
#'@param site.list The list of GRanges objects or GRangesList containing the mapped site coordinates.
#'
#'@return None
#'
#'@import methods
#'@import GenomicRanges
#'
#'@export
#'
check_sites <- function( site.list ){

  if( !is.list( site.list ) & !is( site.list, "GRangesList" ) ){
    stop( "site.list must be a list of GRanges objects or a GRangesList",
          call. = FALSE )
  }

  if( is.null( names( site.list ) ) | any( duplicated( names( site.list ) ) ) ){
    stop( "site.list elements must be uniquely named.",
          call. = FALSE )
  }

  vv <- do.call( c, lapply( site.list, function(x) any( width(x) != 1 ) ) )

  if( any( vv ) ){
    warning( "Not all sites have width = 1.",
             "\n",
             "This can effect quantified overlaps.",
             "\n",
             "In addition, only the start site position is used for sequence logos, meta-feature analysis, etc.",
             "\n",
             "If you think this is a mistake, ensure that only IS datasets were included in your imported data.",
             "\n",
             "Run fix_width() to enforce width = 1.",
             call. = FALSE )
  }
}
