#' Feature Density
#'
#' Calculates the number of features within a defined distance of each site.
#'
#'@param site.list The list or GRangesList containing the mapped site coordinates.
#'@param features The feature-set of interest of interest. Should be a GRanges object.
#'@param genome.obj The genome object of interest.
#'@param win.size The total distance around each site.
#'@param min.overlap The minimum amount of overlap required to call two regions overlapping. Defaults to 1L.
#'@param average Boolean. Whether or not to average the number of features within win.size for each dataset. Defaults to TRUE.
#'@param remove.outliers Whether or not to remove regions that, when expanded, violate genome boundaries.
#'
#'@export
#'
feature_density <- function( site.list, features,  genome.obj, win.size = 1E6, min.overlap = 1L, average = TRUE, remove.outliers = FALSE ){

  fd <- lapply( X = site.list,
                FUN = function(x){
                  expanded <- resize( x, width = win.size, fix = "center", ignore.strand = TRUE )
                  outs <- bound_check( fragments = expanded,
                                       genome.obj = genome.obj,
                                       include.lower = TRUE )
                  if( length(outs) != 0 ){
                    warning( "Expanded coordinates " , paste(outs, collapse=", "),
                             " are out of bounds. To remove these, set \'remove.outliers=TRUE\'.",
                             call. = FALSE )
                    if( isTRUE(remove.outliers) ){
                      expanded <- expanded[-outs]
                      }
                    }

                  ol <- countOverlaps( query = expanded,
                                       subject = features,
                                       minoverlap = min.overlap,
                                       type = "any",
                                       ignore.strand = TRUE )
                  if( isTRUE( average ) ){
                    return( mean( ol ) )
                    } else{
                      return( ol )
                    }
                  }
                )
  if( isTRUE( average ) ){
    return( do.call( c, fd ) )
    } else{
      return( fd )
    }
}
