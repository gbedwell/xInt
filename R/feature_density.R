#' Calculate Feature Density
#'
#' Calculates the number of features within a defined distance of each integration site.
#'
#'@param site.list The list or GRangesList containing the mapped site coordinates.
#'@param features The feature-set of interest. Should be a GRanges object.
#'@param genome.obj The BSgenome object of interest.
#'Used to assess center coordinates and assess outlier regions.
#'@param win.size The total distance around each site.
#'@param current.start The position in the target site duplication currently described by start.
#'This is used for centering the IS coordinates.
#'@param tsd The total length of the target site duplication. This is used for centering the IS coordinates.
#'@param min.overlap The minimum amount of overlap required to call two regions overlapping. Defaults to 1L.
#'@param average Boolean. Whether or not to average the number of features within win.size for each dataset.
#'Defaults to TRUE and returns a single value for each dataset.
#'@param remove.outliers Boolean. Whether or not to remove regions that, when expanded, violate genome boundaries.
#'Defaults to FALSE.
#'
#'@return If average = TRUE, a single vector containing the average feature density value for each dataset.
#'Otherwise, a list of numeric vectors for each dataset containing the number of features within each window.
#'
#'@examples
#'library(BSgenome.Hsapiens.UCSC.hs1)
#'data(sites)
#'data(xobj)
#'feats <- rowRanges(xobj)
#'feature_density(site.list = sites,
#'                features = feats,
#'                genome.obj = Hsapiens)
#'
#'@import methods
#'@import GenomicRanges
#'
#'@export
#'
feature_density <- function( site.list, features,  genome.obj = NULL, win.size = 1E6,
                             current.start = 1, tsd = 5, min.overlap = 1L, average = TRUE,
                             remove.outliers = FALSE ){

  # if( !validObject( site.list ) ){
  #   stop( "site.list is not a valid SiteListObject.",
  #         call. = FALSE )
  # }

  cc <- center_coordinates( site.list = site.list,
                            current = current.start,
                            tsd = tsd,
                            genome.obj = genome.obj )

  fd <- lapply( X = cc,
                FUN = function(x){
                  expanded <- resize( x, width = win.size, fix = "center", ignore.strand = TRUE )

                  if( isTRUE( remove.outliers ) ){

                    if( is.null(genome.obj) ){
                      stop("genome.obj cannot be NULL if remove.outliers is TRUE.")
                    }

                    outs <- bound_check( fragments = expanded,
                                         genome.obj = genome.obj,
                                         include.lower = TRUE )

                    if( length(outs) != 0 ){

                      warning( "Expanded coordinates " , paste(outs, collapse=", "),
                               " contain out of bounds ranges. These coordinates will be removed.",
                               call. = FALSE )

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
    } else{
      names(fd) <- names(site.list)
    }

  return( fd )
}
