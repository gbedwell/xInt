#' Collapse Overlapping Features
#'
#' Collapse overlapping features-of-interest into a single region.
#'
#' @param features A GRanges object holding the features-of-interest.
#' @param ignore.strand Boolean. Whether or not to ignore the strandedness of the features. Defaults to TRUE.
#' @param original.names Boolean. Whether or not keep the original feature names.
#' In the case of overlapping features, the feature names are changed to include them all, separated by hyphens.
#' This can be slow.
#' If FALSE (default), new features are named Feature1 to FeatureN.
#' @param id.col The metadata column containing feature names.
#' Must be provided if original.names = TRUE.
#'
#' @return A GRanges object of collapsed features.
#'
#' @examples
#' data(xobj)
#' feats <- rowRanges(xobj)
#' collapse_features(features = feats)
#'
#' @import GenomicRanges
#'
#' @export
#'
collapse_features <- function(features, ignore.strand = TRUE, original.names = FALSE, id.col = NULL){

  collapsed <- reduce(x = features, ignore.strand = ignore.strand, with.revmap = TRUE)

  if (isTRUE(original.names)){
    if(is.null(id.col) | missing(id.col)){
      stop("ID column must be named to extract feature names.", call. = FALSE)
    }

    cnames <- lapply(X = mcols(collapsed)$revmap,
                     FUN = function(x){
                       paste0(mcols(features)[ , id.col ][ x ], collapse = "-")
                       }
                     )
    mcols(collapsed)$name <- as.character(cnames)
    } else{
      mcols(collapsed)$name <- paste0("Feature", seq_along(collapsed))
    }
  return(collapsed)
}
