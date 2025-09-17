#' Meta-Feature Analysis
#'
#' Calculates the relative position of sites along defined features-of-interest.
#' Returns the relative positions for each site in each dataset.
#'
#' @param sites The list or GRangesList containing IS coordinates.
#' @param features The features-of-interest. Must be a GRanges object.
#' @param collapse Boolean. Whether or not to collapse the output into a single data frame.
#' Defaults to TRUE. For compatibility with plot_metafeature(), this must be TRUE.
#' @param conditions The condition for each dataset. Must be same length as sites.
#' 
#' @return A data frame or a list of data frames.
#'
#' @examples
#' data(sites2)
#' data(xobj)
#' feats <- rowRanges(xobj)
#'
#' metafeature(sites = sites2,
#'             features = feats)
#'
#' @import GenomicRanges
#' @import S4Vectors
#'
#' @export
#'
metafeature <- function(sites,
                        features,
                        collapse = TRUE,
                        conditions = NULL){

  if(!validObject(sites)){
    stop("sites is not a valid SiteList object.",
         call. = FALSE)
  }
  
  if(!is.null(conditions) && length(conditions) != length(sites)){
    stop("conditions must have the same length as sites.",
         call. = FALSE)
  }

  pos.ll <- lapply(
    X = sites,
    FUN = function(x){
      gr <- x
      ol <- findOverlaps(
        query = gr,
        subject = features,
        minoverlap = 1L,
        type = "any",
        ignore.strand = TRUE
        )

      ol.sites <- gr[queryHits(ol)]
      ol.feats <- features[subjectHits(ol)]

      df <- data.frame(
        chr = as.character(seqnames(ol.sites)),
        site = start(ranges(ol.sites)),
        site.strand = as.character(strand(ol.sites)),
        feature.start = start(ranges(ol.feats)),
        feature.end = end(ranges(ol.feats)),
        feature.strand = as.character(strand(ol.feats))
        )

      df$rel.position <- (df$site - df$feature.start) / (df$feature.end - df$feature.start)
      df[df$feature.strand == "-", ]$rel.position <- 1 - df[df$feature.strand == "-", ]$rel.position

      return(df)
      }
    )

  # Add dataset names and condition information
  pos.ll <- Map(function(x, y) { 
    x$dataset <- y
    if(!is.null(conditions)) {
      dataset.to.condition <- setNames(conditions, names(sites))
      x$condition <- dataset.to.condition[y]
    } else {
      x$condition <- x$dataset
    }
    return(x)
  }, pos.ll, names(sites))
  
  # Always combine by condition if conditions are provided
  if(!is.null(conditions)) {
    if(isTRUE(collapse)) {
      pos.ll <- do.call(rbind, pos.ll)
      rownames(pos.ll) <- NULL
    } else {
      all.datasets <- do.call(rbind, pos.ll)
      pos.ll <- split(all.datasets, all.datasets$condition)
    }
  } else if(isTRUE(collapse)) {
    pos.ll <- do.call(rbind, pos.ll)
    rownames(pos.ll) <- NULL
  }
  
  return(pos.ll)
}