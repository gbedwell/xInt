#' Calculate Feature Density
#'
#' Calculates the number of features within a defined distance of each integration site.
#'
#'@param site.list The list or GRangesList containing the mapped site coordinates.
#'@param features The feature-set of interest. Should be a GRanges object.
#'@param genome.obj The genome object of interest.
#'@param win.size The total distance around each site.
#'@param min.overlap The minimum amount of overlap required to call two regions overlapping. Defaults to 1L.
#'@param average Boolean. Whether or not to average the number of features within win.size for each dataset. Defaults to TRUE and returns a single value for each dataset.
#'@param remove.outliers Boolean. Whether or not to remove regions that, when expanded, violate genome boundaries. Defaults to FALSE.
#'
#'@return If average = TRUE, a single vector containing the average feature density value for each dataset. Otherwise, a list of numeric vectors for each dataset containing the number of features within each window.
#'
#'@import methods
#'@import GenomicRanges
#'
#'@export
#'
feature_density <- function( site.list, features,  genome.obj, win.size = 1E6, min.overlap = 1L,
                             average = TRUE, remove.outliers = FALSE ){

  check_sites( site.list )

  fd <- lapply( X = site.list,
                FUN = function(x){
                  expanded <- resize( x, width = win.size, fix = "center", ignore.strand = TRUE )

                  if( remove.outliers ){
                    outs <- bound_check( fragments = expanded,
                                         genome.obj = genome.obj,
                                         include.lower = TRUE )
                    if( length(outs) != 0 ){
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
    fd <- do.call( c, fd )
    names(fd) <- names(site.list)
    return( fd )
    } else{
      names(fd) <- names(site.list)
      return( fd )
    }
}
