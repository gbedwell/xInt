#' Make xIntNearest Object
#'
#' Creates an object that stores the distance to and identity of the nearest feature 
#' for each integration site across multiple samples.
#'
#' @param sites A SiteList object.
#' @param features A GRanges object for the feature-set of interest.
#' @param conditions The infection condition for each dataset.
#' Must be same length as site.list.
#' @param condition.levels The factor levels of conditions.
#' All entries should be present in conditions.
#' Length should be the number of unique conditions.
#' @param id.col The ID column in the features GRanges object.
#' Must be provided.
#' Default "name".
#' @param ignore.strand Logical indicating whether to ignore strand information when calculating distances.
#' Default is TRUE.
#' @param ... Used for adding additional columns to the object metadata.
#'  Non-numeric and non-double vectors are automatically converted to a factor.
#' Numeric and double vectors are left as-is.
#' If numeric or double vectors are intended to be factors, enter them wrapped in factor().
#' Intended for batch information, etc.
#'
#' @return An xIntNearest object containing:
#'  \itemize{
#'    \item sites: A named list of data frames, one per sample, with site coordinates and nearest feature information
#'    \item summary: A named list of summary statistics for each sample
#'    \item summary_df: A data frame with summary statistics for all samples
#'    \item features: The original features GRanges object
#'    \item metadata: Additional metadata about the analysis
#'  }
#'
#' @examples
#' data(sites)
#' sites <- sites[!names(sites) %in% "C1"]
#' data(xobj)
#' feats <- rowRanges(xobj)
#' make_xIntNearest(sites = sites,
#'                 features = feats,
#'                 conditions = c(rep("A",4),rep("B",5)),
#'                 condition.levels = c("A","B"))
#'
#' @import GenomicRanges
#' @import S4Vectors
#'
#' @export
#'
make_xIntNearest <- function(sites,
                             features,
                             conditions,
                             condition.levels,
                             id.col = "name",
                             ignore.strand = TRUE,
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

  site.list <- as.list(sites)
  sample.names <- names(site.list)

  if(length(sample.names) != length(conditions)){
    stop("The number of sample names does not match the number of conditions.
         There must be 1 condition given for each dataset in site.list.",
         call. = FALSE)
  }

  if(!all(conditions %in% condition.levels)){
    stop("Not all conditions found in condition levels.",
         call. = FALSE)
  }

  # Calculate distances and nearest features
  nearest.results <- lapply(
    X = site.list,
    FUN = function(x){
      # Find nearest feature for each site
      nearest.idx <- nearest(x, features, ignore.strand = ignore.strand)
      
      # Get distances
      distances <- distanceToNearest(x, features, ignore.strand = ignore.strand)
      
      # Get feature IDs
      feature.ids <- mcols(features)[nearest.idx, id.col]
      feature.starts <- start(features)[nearest.idx]
      feature.ends <- end(features)[nearest.idx]
      feature.strands <- as.character(strand(features))[nearest.idx]
      
      # Create data frame with site info, distance, and nearest feature ID
      site.df <- data.frame(
        seqnames = as.character(seqnames(x)),
        start = start(x),
        end = end(x),
        strand = as.character(strand(x)),
        distance = mcols(distances)$distance,
        nearest.feature = feature.ids,
        feature.start = feature.starts,
        feature.end = feature.ends,
        feature.strand = feature.strands
      )
      
      return(site.df)
    }
  )
  
  # Preserve the names from site.list
  names(nearest.results) <- sample.names
  
  # Add sample and condition information to each result
  for (i in seq_along(nearest.results)) {
    nearest.results[[i]]$sample <- sample.names[i]
    nearest.results[[i]]$condition <- conditions[i]
  }
  
  # Create summary statistics for each sample
  sample.stats <- lapply(
    X = nearest.results,
    FUN = function(x){
      data.frame(
        total.sites = nrow(x),
        mean.distance = mean(x$distance),
        median.distance = median(x$distance),
        min.distance = min(x$distance),
        max.distance = max(x$distance)
      )
    }
  )
  
  # Preserve the names from site.list
  names(sample.stats) <- sample.names
  
  # Create summary data frame
  summary.df <- do.call(rbind, sample.stats)
  summary.df <- data.frame(
    sample = sample.names,
    summary.df,
    condition = factor(conditions, levels = condition.levels)
  )
  rownames(summary.df) <- sample.names
  
  # Add additional columns
  other.cols <- list(...)
  
  for(vec in names(other.cols)){
    if(is.numeric(other.cols[[vec]]) ||
       is.double(other.cols[[vec]]) ||
       is.factor(other.cols[[vec]])){
      summary_df[[vec]] <- other.cols[[vec]]
    } else{
      summary_df[[vec]] <- factor(other.cols[[vec]])
    }
  }
  
  # Create the xIntNearest object as an S3 class
  result <- list(
    sites = nearest.results,
    summary.df = summary.df
  )
  
  class(result) <- "xIntNearest"

  validation <- validate_xIntNearest(result)
  if (is.character(validation)) {
    stop(validation, call. = FALSE)
  }
  
  return(result)
}
