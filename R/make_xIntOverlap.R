#' Make xIntOverlap Object
#'
#' Creates a ranged summarized experiment object compatible with other xInt functions.
#' Stores summary and global information for each sample in column data.
#' Range information is stored in row data.
#'
#'@param sites A list of GRanges objects or a GRangesList holding IS coordinates.
#'@param features A GRanges object for the feature-set of interest.
#'@param conditions The condition for each dataset. Must be same length as sites.
#'@param condition.levels The factor levels of conditions.
#'All entries should be present in conditions.
#'Length should be the number of unique conditions.
#'@param min.overlap The minimum amount of overlap to count two regions as overlapping.
#'Defaults to 1.
#'@param id.col The ID column in the features GRanges object.
#'Must be provided.
#'Default "name".
#' @param expand Numeric vector of length 2 specifying the number of base pairs to expand the given ranges up/downstream.
#' Defaults to c(0,0).
#'@param ... Used for adding additional columns to the object column data.
#'Non-numeric and non-double vectors are automatically converted to a factor.
#'Numeric and double vectors are left as-is.
#'If numeric or double vectors are intended to be factors, enter them wrapped in factor().
#'Intended for batch information, etc.
#'
#'@return A SummarizedExperiment object containing local and global feature overlap
#'information for each dataset in sites.
#'
#'@examples
#'data(sites)
#'sites <- sites[!names(sites) %in% "C1"]
#'data(xobj)
#'feats <- rowRanges(xobj)
#'make_xIntOverlap(sites = sites,
#'                features = feats,
#'                conditions = c(rep("A",4),rep("B",5)),
#'                condition.levels = c("A","B"))
#'
#'@import GenomicRanges
#'@import SummarizedExperiment
#'@import S4Vectors
#'@importFrom tidyr pivot_wider
#'@importFrom tibble column_to_rownames
#'
#'@export
#'
make_xIntOverlap <- function(sites,
                             features,
                             conditions,
                             condition.levels,
                             min.overlap = 1L,
                             id.col = "name",
                             expand = c(0,0),
                             ... ){

  if(!id.col %in% names(mcols(features))){
    stop("features must have an identifier column holding identifiers for each annotated region.",
         "\n",
         "The name of this column should be provided in the id.col argument.",
         call. = FALSE)
  }

  if(!validObject(sites)){
    stop("sites is not a valid SiteList object.",
         call. = FALSE )
  }

  # sites <- as.list(sites)
  sample.names <- names(sites)

  if(length(sample.names) != length(conditions)){
    stop("The number of sample names does not match the number of conditions.
         There must be 1 condition given for each dataset in sites.",
         call. = FALSE)
  }

  if(!all(conditions %in% condition.levels)){
    stop("Not all conditions found in condition levels.",
         call. = FALSE)
  }

  if(any(abs(expand) > 0)) {
    features <- expand_features(
      features = features,
      upstream = expand[1], 
      downstream = expand[2]
    )
  }

  hits <- lapply(
    X = sites,
    FUN = function(x){
    findOverlaps(
      query = features,
      subject = x,
      minoverlap = min.overlap,
      type = "any",
      ignore.strand = TRUE
      )
    }
  )

  frac.overlap <- lapply(
    X = hits,
    FUN = function(x){
      total.sites <- length(countSubjectHits(x))
        overlapping <- length(unique(subjectHits(x)))
        frac.overlap <- overlapping/total.sites
        multioverlap <- length(countSubjectHits(x)[countSubjectHits(x) > 1])

        if(multioverlap != 0){
          ol.check <- length(
            countSubjectHits(x)[countSubjectHits(x) != 0]) == overlapping

          if(!isTRUE(ol.check)){
            stop("Error in counting feature overlaps.", call. = FALSE)
            }
          }

        data.frame(
          total.sites = total.sites,
          overlapping.sites = overlapping,
          fraction.overlap = frac.overlap
          )
        }
      )

  frac.overlap <- do.call(rbind, frac.overlap)

  frac.overlap <- data.frame( 
    sample = sample.names,
    frac.overlap,
    condition = factor(conditions, levels = condition.levels
    )
  )

  other.cols <- list(...)

  for(vec in names(other.cols)){
    if(is.numeric( other.cols[[vec]]) ||
       is.double( other.cols[[vec]]) ||
       is.factor( other.cols[[vec]])){
      frac.overlap[[vec]] <- other.cols[[vec]]
    } else{
      frac.overlap[[vec]] <- factor(other.cols[[vec]])
    }
  }

  feature.counts <- lapply(
    X = hits,
    FUN = function(x){
      counts <- c(countQueryHits(x))
      ids <- c(mcols(features)[, names(mcols(features)) %in% id.col])
      ids <- make.unique(ids, sep = ".dup")
      return(data.frame(ids, counts))
      }
    )

  feature.counts <- Map(cbind, feature.counts, sample = sample.names)

  # TO-DO: Rewrite this to remove tidyverse functions.
  count.mat <- as.matrix(
    do.call( rbind, c( feature.counts, make.row.names = FALSE ) ) |>
      tidyr::pivot_wider( names_from = sample, values_from = counts ) |>
      tibble::column_to_rownames( var = "ids" )
    )
  
  count.mat[is.na(count.mat)] <- 0L

  # Ensure count values are integers
  storage.mode(count.mat) <- "integer"

  xint.obj <- SummarizedExperiment( 
    assays = list(counts = count.mat),
    colData = frac.overlap,
    rowRanges = features 
    )

  xint.obj <- new("xIntOverlap", xint.obj)

  return(xint.obj)
}
