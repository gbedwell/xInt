#' Combine Column Data
#'
#' Combines colData from multiple xIntObjects into a single data frame.
#' A new 'annotation' column is added to the resulting data frame,
#' holding the name of the annotations used to create each xIntObject.
#'
#'@param ... The xIntObjects of interest.
#'@param annotations A character vector containing the names of the feature-sets used to create each xIntObject.
#'
#'@return A data frame containing the colData from each xIntObject.
#'
#'@examples
#'data(xobj)
#'yobj <- xobj
#'join_colData(xobj, yobj, annotations = c("Anno1", "Anno2"))
#'
#'@import SummarizedExperiment
#'
#'@export
#'
join_colData <- function( ..., annotations ){

  if( !all( do.call( c, lapply( list( ... ), validObject ) ) ) ){
    stop( "All provided objects must be valid xIntObjects.",
          call. = FALSE )
  }

  if( length( list( ... ) ) != length( annotations ) ){
    stop( "The number of annotations must equal the number of provided objects.",
          call. = FALSE )
  }

  ll <- lapply(
    X = list( ... ),
    FUN = function(x){
      dat <- colData(x)
      rownames(dat) <- NULL
      return( dat )
    }
  )

  ll <- Map( function(x, y) { cbind( x, annotation = y ) }, ll, annotations )

  df <- data.frame( do.call( rbind, ll ) )

  return( df )
}
