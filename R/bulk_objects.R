#' Make xInt Objects in Batch
#'
#' Make xIntOverlap or xIntNearest objects in bulk from a list of feature-sets.
#'
#' @param site A SiteList object.
#' @param conditions The condition for each dataset. Must be same length as sites.
#' @param condition.levels The factor levels of conditions.
#' All entries should be present in conditions.
#' Length should be the number of unique conditions.
#' @param min.overlap The minimum amount of overlap to count two regions as overlapping.
#' Defaults to 1.
#' @param id.col The ID column in the features GRanges object. Default "name".
#' @param features The list containing to different feature-sets.
#' @param type The type of xInt object to construct. 
#' Can be "overlap" for xIntOverlap (default) or "nearest" for xIntNearest.
#' @param feature.names A character vector of names for each feature in features.
#' When NULL (default), features must be named.
#' @param ric.dat A separate SiteList object for random integration control sites.
#' Can be NULL (default).
#'
#' @return A list of xInt objects.
#'
#' @export

bulk_objects <- function(sites, conditions, condition.levels, features,
                         type = c("overlap", "nearest"), feature.names = NULL,
                         min.overlap = 1, id.col = "name", ric.dat = NULL){

  if(!validObject(sites)){
    stop("sites is not a valid SiteList object.",
         call. = FALSE )
  }

  if(!is.list(features)) {
    stop("features must be a list.", call. = FALSE)
  }

  type = match.arg(type)
  
  if(is.null(feature.names)){
    if(is.null(names(features))){
      stop("features must be named if feature.names is NULL.")
    }
    else{
      feature.names <- names(features)
    }
  }

  if(type == "overlap") {
    overlap.objects <- list()

    for (i in seq_along(feature.names)) {
      overlap.objects[[i]] <- make_xIntOverlap(
        sites = sites,
        features = features[[i]],
        conditions = conditions,
        condition.levels = condition.levels,
        min.overlap = min.overlap,
        id.col = id.col
      )
    }
    
    names(overlap.objects) <- feature.names
    
    dat <- do.call(join_colData, c(overlap.objects, list(annotations = feature.names)))
    
    if(!is.null(ric.dat)){
      dat <- rbind(dat, ric.dat)
      dat$condition <- factor(dat$condition, levels = c(condition.levels, "RIC"))
    } else{
      dat$condition <- factor(dat$condition, levels = condition.levels)
    }
    
    dat$annotation <- factor(dat$annotation, levels = feature.names)
    
    result <- c(overlap.objects, list(joined = dat))
    } else if(type == "nearest") {
      nearest.objects <- list()

      for (i in seq_along(feature.names)) {
        nearest.objects[[i]] <- make_xIntNearest(
          sites = sites,
          features = features[[i]],
          conditions = conditions,
          condition.levels = condition.levels
        )
      }

      names(nearest.objects) <- feature.names

      result <- nearest.objects
  }
  
  return(result)
}