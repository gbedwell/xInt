#' xIntOverlap-class
#' @examples
#' data(xobj)
#' validObject(xobj)
#'
#' @import methods
#' @export
setClass("xIntOverlap",
         contains = "RangedSummarizedExperiment")

setValidity(
  "xIntOverlap",
  function(object){
    if(!("counts" %in% assayNames(object))){
      return("The assay slot must contain a matrix named 'counts'.")
    }

    if(any(is.na(assays(object)))){
      return("Count data cannot contain NA.")
    }

    if(!is.integer(assays(object)$counts)){
      return("Count data must be integers.")
    }

    if(any(assays(object)$counts < 0)){
      return("Count data must be positive.")
    }

    if(!all(c("sample", "total.sites", "overlapping.sites", "fraction.overlap", "condition") %in%
              names(colData(object)))){

      msg <- paste("colData must contain",
                    paste(c("sample", "total.sites", "overlapping.sites", "fraction.overlap", "condition"),
                           collapse = ", "),
                    collapse = " ")
      return(msg)
    }

    if(!is.integer(colData(object)$total.sites)){
      return("'total.sites' must be integers.")
    }

    if(!is.integer(colData(object)$overlapping.sites)){
      return("'overlapping.sites' must be integers.")
    }

    if(!is.numeric(colData(object)$fraction.overlap)){
      return("'fraction.overlap' must be numeric.")
    }

    if(!is.factor(colData(object)$condition)){
      return("'condition' column must be a factor.")
    }

    if(!all(colData(object)$sample == colnames(object))){
      return("count matrix column names and 'sample' column in colData must match.")
    }

    return(TRUE)

  }
)

#' Constructor for xIntOverlap
#'
#' @param counts A matrix containing integer count data with row names for features and column names for samples.
#' @param colData A DataFrame containing the sample metadata with columns: 'sample', 'total.sites',
#' 'overlapping.sites','fraction.overlap', and 'condition'.
#' @param ... Additional arguments passed to the RangedSummarizedExperiment constructor.
#'
#' @return An xIntOverlap object.
#' @export
xIntOverlap <- function(counts, colData, ...) {
  # Check that counts is a matrix and contains integers
  if (!is.matrix(counts) || !is.integer(counts)) {
    stop("counts must be an integer matrix.")
  }

  # Check that colData is a DataFrame and has the necessary columns
  if (!is.data.frame(colData) || !all(c("sample", "total.sites", "overlapping.sites", "fraction.overlap", "condition") %in% names(colData))) {
    stop("colData must be a DataFrame containing the required columns: sample, total.sites, overlapping.sites, fraction.overlap, condition.")
  }

  # Ensure the sample names in colData match the column names in counts
  if (!all(colData$sample == colnames(counts))) {
    stop("Sample names in colData must match the column names in counts.")
  }

  # Create a RangedSummarizedExperiment object
  rse <- RangedSummarizedExperiment(assays = list(counts = counts), colData = colData, ...)

  # Create and return a xIntOverlap object
  new("xIntOverlap", rse)
}

#' SiteList-class
#'
#' @examples
#' sites <- SiteList(list(GRanges(seqnames="chr1", ranges=IRanges(start=45, end=45))))
#' validObject(sites)
#'
#' @import methods
#' @import GenomicRanges
#' @export
setClass(
  "SiteList",
  slots = list(
    sites = "ANY"
    )
  )

setValidity(
  "SiteList",
  function(object) {
    if (!is(object@sites, "GRangesList") &&
        !is.list(object@sites)) {
      return("sites slot must be a GRangesList or a list of GRanges objects.")
    }

    if (is.list(object@sites)) {
      if (!all(sapply(object@sites, function(x) is(x, "GRanges")))) {
        return("All elements of the list must be GRanges objects.")
      }
      
      # Check that all ranges have width 1
      widths.check <- all(sapply(object@sites, function(gr) {
        all(width(gr) == 1)
      }))
      
      if (!widths.check) {
        return("All ranges in all GRanges objects must have a width of 1.")
      }
    } else if (is(object@sites, "GRangesList")) {
      # Check that all ranges in GRangesList have width 1
      if (!all(unlist(width(object@sites)) == 1)) {
        return("All ranges in the GRangesList must have a width of 1.")
      }
    }
    
    return(TRUE)
  }
)

#' Constructor for SiteList
#'
#' @param sites A GRangesList or a list of GRanges objects
#' @return A SiteList object
#' @examples
#' SiteList(list(GRanges(seqnames="chr1", ranges=IRanges(start=45, end=45))))
#'
#' @export
SiteList <- function(sites) {
  if(!is(sites, "GRangesList") && !is.list(sites)){
    stop("sites must be a GRangesList or a list of GRanges objects.")
  }
  
  if(is.list(sites) && !all(sapply(sites, function(x) is(x, "GRanges")))){
    stop("All elements of the list must be GRanges objects.")
  }
  
  # Check that all ranges have width 1
  if(is.list(sites)) {
    widths_check <- all(sapply(sites, function(gr) {
      all(width(gr) == 1)
    }))
    
    if(!widths_check) {
      stop("All ranges in all GRanges objects must have a width of 1.")
    }
  } else if(is(sites, "GRangesList")) {
    if(!all(unlist(width(sites)) == 1)) {
      stop("All ranges in the GRangesList must have a width of 1.")
    }
  }
  
  new("SiteList", sites = sites)
}

#' @export
setMethod("[", "SiteList", function(x, i, j, ..., drop = TRUE) {
  if (missing(j)) {
    result <- x@sites[i]
    if (is.list(result) && length(result) > 0) {
      return(SiteList(result))
    } else {
      return(result)
    }
  } else {
    stop("Invalid subsetting operation for SiteList")
  }
})

#' @export
setMethod("[[", "SiteList", function(x, i, j, ...) {
  if (missing(j)) {
    return(x@sites[[i]])
  } else {
    stop("Invalid subsetting operation for SiteList")
  }
})

#' @export
setMethod("$", "SiteList", function(x, name) {
  return(x@sites[[name]])
})

#' @export
setMethod("length", "SiteList", function(x) {
  return(length(x@sites))
})

#' @export
setMethod("names", "SiteList", function(x) {
  return(names(x@sites))
})

#' @export
setReplaceMethod("names", "SiteList", function(x, value) {
  names(x@sites) <- value
  return(x)
})

#' @export
setMethod("show", "SiteList", function(object) {
  cat("SiteList object with", length(object), "elements\n")
  if (length(object) > 0) {
    if (!is.null(names(object))) {
      cat("Names:", paste(names(object), collapse=", "), "\n")
    }
    cat("Each element is a GRanges object\n")
    for (i in 1:min(3, length(object))) {
      cat("Element", i, ":", length(object[[i]]), "sites\n")
    }
    if (length(object) > 3) {
      cat("...\n")
    }
  }
})

#' @export
setAs("SiteList", "list", function(from) {
  return(from@sites)
})

#' @export
setMethod("as.list", "SiteList", function(x, ...) {
  return(x@sites)
})

#' @export
setMethod("lapply", "SiteList", function(X, FUN, ...) {
  result <- lapply(X@sites, FUN, ...)
  if (all(sapply(result, function(x) is(x, "GRanges")))) {
    return(SiteList(result))
  } else {
    return(result)
  }
})

#' @export
setMethod("sapply", "SiteList", function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  return(sapply(X@sites, FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES))
})

#' xIntNearest class
#'
#' @description
#' An S3 class that stores the distance to and identity of the nearest feature 
#' for each integration site across multiple samples.
#'
#' @details
#' The xIntNearest object contains two main components:
#' \itemize{
#'   \item sites: A named list of data frames, one per sample, with site coordinates and nearest feature information
#'   \item summary.df: A data frame with summary statistics for all samples
#' }
#'
#' @param sites A named list of data frames, one per sample, with site coordinates and nearest feature information
#' @param summary.df A data frame with summary statistics for all samples
#'
#' @return An object of class \code{xIntNearest}
#'
#' @examples
#' data(xnear)
#'
#' @export
xIntNearest <- function(sites, summary.df) {
  # Validate inputs
  if (!is.list(sites) || length(sites) == 0) {
    stop("'sites' must be a non-empty list of data frames.")
  }
  
  if (!all(sapply(sites, is.data.frame))) {
    stop("All elements in 'sites' must be data frames.")
  }
  
  if (!is.data.frame(summary.df)) {
    stop("'summary.df' must be a data frame.")
  }
  
  # Check required columns in summary.df
  required_cols <- c("sample", "total.sites", "mean.distance", "median.distance", 
                     "min.distance", "max.distance", "condition")
  if (!all(required_cols %in% colnames(summary.df))) {
    stop("'summary.df' must contain the following columns: ", 
         paste(required_cols, collapse = ", "))
  }
  
  # Create the xIntNearest object
  result <- list(
    sites = sites,
    summary.df = summary.df
  )
  
  class(result) <- "xIntNearest"
  
  return(result)
}

validate_xIntNearest <- function(object) {
  if (!is.list(object)) return("Object must be a list")
  if (!all(c("sites", "summary.df") %in% names(object))) 
    return("Object must contain 'sites' and 'summary.df' elements")
  # Add more validation as needed
  return(TRUE)
}

#' @export
print.xIntNearest <- function(x, ...) {
  class(x) <- NULL
  print(x, ...)
  invisible(x)
}