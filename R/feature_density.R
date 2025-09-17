#' Calculate Feature Density
#'
#' Calculates the number of defined features within a defined distance of each integration site.
#'
#' @param sites A SiteList object.
#' @param features The feature-set of interest. Should be a GRanges object.
#' @param genome.obj The BSgenome object of interest.
#' Used to assess center coordinates and assess outlier regions.
#' @param win.size The total distance around each site.
#' @param current.start The position in the target site duplication currently described by start.
#' This is used for centering the IS coordinates.
#' @param tsd The total length of the target site duplication. This is used for centering the IS coordinates.
#' @param min.overlap The minimum amount of overlap required to call two regions overlapping. Defaults to 1L.
#' @param average Boolean. Whether or not to average the number of features within win.size for each dataset.
#' Defaults to TRUE and returns a single value for each dataset.
#' @param remove.outliers Boolean. Whether or not to remove regions that, when expanded, violate genome boundaries.
#' Defaults to FALSE.
#' @param conditions The condition for each dataset. Must be the same length and in the same order as sites.
#'
#' @return If average = TRUE, a single vector containing the average feature density value for each dataset.
#' Otherwise, a list of numeric vectors for each dataset containing the number of features within each window.
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hs1)
#' data(sites)
#' data(xobj)
#' feats <- rowRanges(xobj)
#' feature_density(sites = sites,
#'                 features = feats,
#'                 genome.obj = Hsapiens)
#'
#' @import methods
#' @import GenomicRanges
#'
#' @export
#'
feature_density <- function(sites, features,  genome.obj = NULL, win.size = 1E6,
                            current.start = 1, tsd = 5, min.overlap = 1L, average = TRUE,
                            remove.outliers = FALSE, conditions = NULL){

  if(!validObject(sites)){
    stop("sites is not a valid SiteList object.",
         call. = FALSE )
  }

  cc <- center_coordinates(sites = sites,
                           current = current.start,
                           tsd = tsd,
                           genome.obj = genome.obj)

  fd <- lapply( 
    X = cc,
    FUN = function(x){
      expanded <- resize(x, width = win.size, fix = "center", ignore.strand = TRUE)

      if(isTRUE(remove.outliers)){

        if(is.null(genome.obj)){
          stop("genome.obj cannot be NULL if remove.outliers is TRUE.")
        }

        outs <- bound_check(fragments = expanded,
                            genome.obj = genome.obj,
                            include.lower = TRUE )

        if(length(outs) != 0){

          warning("Expanded coordinates " , paste(outs, collapse=", "),
                  " contain out of bounds ranges. These coordinates will be removed.",
                  call. = FALSE )

          expanded <- expanded[-outs]
          }
        }

      ol <- countOverlaps(query = expanded,
                          subject = features,
                          minoverlap = min.overlap,
                          type = "any",
                          ignore.strand = TRUE )
      if(isTRUE(average)){
        return(mean(ol))
        } else{
          return(ol)
        }
      }
    )

  if(isTRUE(average)){
    fd <- do.call(c, fd)
    names(fd) <- names(sites)

    if(!is.null(conditions)) {
      if(length(conditions) != length(fd)) {
        stop("Length of conditions must match the number of samples in sites.",
             call. = FALSE)
      }

      result.df <- data.frame(
        sample = names(fd),
        condition = conditions,
        avg.density = as.numeric(fd)
      )
      
      return(result.df)

    }
    } else{
      names(fd) <- names(sites)
    }

  return(fd)
}
