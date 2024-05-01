#' Collapse Overlapping Features
#'
#' Collapse features-of-interest to eliminate overlapping features that might result in multiple counting of individual integration sites.
#'
#'@param features A GRanges object holding the features-of-interest.
#'@param ignore.strand Logical. Whether or not to ignore the strandedness of the features. Defaults to TRUE.
#'@param original.names Logical. Whether or not keep the original feature names. In the case of overlapping features, the feature names are changed to include them all, separated by hyphens. This can be slow. If FALSE, new features are named Feature 1 to Feature N.
#'@param id.col The metadata column containing feature names. Must be provided if original.names = TRUE. Can be NULL.
#'
#'@return A GRanges object of collapsed features.
#'
#'@import GenomicRanges
#'
#'@export
#'
collapse_features <- function( features, ignore.strand = TRUE, original.names = FALSE, id.col = NULL ){

  collapsed <- reduce( x = features, ignore.strand = ignore.strand, with.revmap = TRUE )

  if ( isTRUE( original.names ) ){
    if( is.null( id.col ) ){
      stop( "ID column must be named to extract feature names.", call. = FALSE )
    }

    cnames <- lapply( X = mcols( collapsed )$revmap,
                      FUN = function(x){
                        paste0( mcols( features )[ , id.col ][ x ], collapse = "-" )
                      } )
    mcols( collapsed )$name <- as.character( cnames )
  } else{
    mcols( collapsed )$name <- paste0( "Feature", seq_along( collapsed ) )
  }
  return( collapsed )
}
