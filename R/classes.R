#'xIntOverlap-class
#'@examples
#'data(xobj)
#'validObject(xobj)
#'
#'@import methods
#'@export
setClass( "xIntOverlap",
          contains = "RangedSummarizedExperiment" )

setValidity(
  "xIntOverlap",
  function( object ){

    if ( !( "counts" %in% assayNames( object ) ) ){
      return( "The assay slot must contain a matrix named 'counts'." )
    }

    if ( any( is.na( assays( object ) ) ) ){
      return( "Count data cannot contain NA." )
    }

    if ( !is.integer( assays( object )$counts ) ){
      return( "Count data must be integers." )
    }

    if ( any( assays( object )$counts < 0 ) ){
      return( "Count data must be positive." )
    }

    if( !all( c( "sample", "total.sites", "overlapping.sites", "fraction.overlap", "condition" ) %in%
              names( colData( object ) ) ) ){

      msg <- paste( "colData must contain",
                    paste( c( "sample", "total.sites", "overlapping.sites", "fraction.overlap", "condition" ),
                           collapse = ", " ),
                    collapse = " " )
      return( msg )
    }

    if( !is.integer( colData( object )$total.sites ) ){
      return( "'total.sites' must be integers." )
    }

    if( !is.integer( colData( object )$overlapping.sites ) ){
      return( "'overlapping.sites' must be integers." )
    }

    if( !is.numeric( colData( object )$fraction.overlap ) ){
      return( "'fraction.overlap' must be numeric." )
    }

    if( !is.factor( colData( object )$condition ) ){
      return( "'condition' column must be a factor.")
    }

    if( !all( colData( object )$sample == colnames( object ) ) ){
      return( "count matrix column names and 'sample' column in colData must match." )
    }

    return( TRUE )

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

  # Create and return a TestOverlap object
  new("TestOverlap", rse)
}


#' SiteList-class
#'
#' @examples
#' sites <- SiteList(list(GRanges(seqnames="chr1", ranges=IRanges(start=1, end=100))))
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
    }

    return(TRUE)
  }
)

SiteList <- function(sites) {
  if (!is(sites, "GRangesList") && !is.list(sites)) {
    stop("sites must be a GRangesList or a list of GRanges objects.")
  }
  new("SiteList", sites = sites)
}




