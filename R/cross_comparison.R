#' Calculate Cross-Sample Overlap
#'
#' Calculates the number of overlapping sites between each sample in a SiteListObject.
#'
#'@param site.list A SiteListObject object.
#'@param unique.sites Boolean. Whether or not to enforce unique sites before calculating overlap.
#'Defaults to TRUE.
#'@param overlap.coef Boolean. Whether or not to return the Szymkiewicz–Simpson overlap coefficient.
#'Defaults to TRUE.
#'
#'@return A matrix of counts denoting the number of overlapping sites between each sample.
#'
#'@examples
#'cross_comparison(site.list = sites)
#'
#'@import GenomicRanges
#'
#'@export
#'
cross_comparison = function( site.list, unique.sites = TRUE, overlap.coef = TRUE ){

  if( !validObject( site.list ) ){
    stop( "site.list is not a valid SiteListObject.",
          call. = FALSE )
  }

  site.count <- lapply(
    X = seq_along( site.list ),
    FUN = function(x){
      tmp.counts <- lapply(
        X = seq_along( site.list ),
        FUN = function(y){
          a <- site.list[[y]]
          b <- site.list[[x]]

          if( isTRUE( unique.sites ) ){
            res <- sum( unique(a) %in% unique(b) )
            } else{
            res <- sum( a %in% b )
            }

          if( isTRUE( overlap.coef ) ){
            if( isTRUE( unique.sites ) ){
              res <- res / min(length(unique(a)), length(unique(b)))
            } else{
              res <- res / min(length(a), length(b))
            }

          }
          return(res)
        }
      )
      do.call(c, tmp.counts)
    }
  )
  mat <- do.call(rbind, site.count)
  colnames(mat) <- names( site.list )
  rownames(mat) <- names( site.list )
  return(mat)
}
