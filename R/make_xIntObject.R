#' Make xIntObject
#'
#' Creates a ranged summarized experiment object compatible with other xInt functions. Stores summary information for each sample in column data. Range information is stored in row data.
#'
#'@param site.list The list of integration site datasets
#'@param features The feature-set of interest.
#'@param conditions The infection condition for each dataset. Must be same length as site.list.
#'@param condition.levels The factor levels of conditions. Length should be the number of unique conditions.
#'@param min.overlap The minimum amount of overlap to count two regions as overlapping
#'@param id.col The ID column in the features GRanges object. Default "name".
#'
#'@return A SummarizedExperiment object containing local and global feature overlap information for each dataset in site.list.
#'
#'@import GenomicRanges
#'@import SummarizedExperiment
#'@import S4Vectors
#'@importFrom tidyr pivot_wider
#'@importFrom tibble column_to_rownames
#'
#'@export
#'
make_xIntObject <- function( site.list,
                         features,
                         conditions,
                         condition.levels,
                         min.overlap = 1L,
                         id.col = "name" ){

  # TO-DO: Add a check for feature rownames. They must either be NULL or match the output count matrix rownames.

  if( !id.col %in% names( mcols( features ) ) ){
    stop( "features must have an identifier column holding identifiers for each annotated region.",
          "\n",
          "The name of this column should be provided in the id.col argument.",
          call. = FALSE )
  }

  check_sites( site.list )

  sample.names <- names( site.list )

  if( length( sample.names ) != length( conditions ) ){
    stop( "The number of sample names does not match the number of conditions. There must be 1 condition given for each dataset in site.list.",
          call. = FALSE )
  }

  if( !all( conditions %in% condition.levels ) ){
    stop( "Not all conditions found in condition levels.",
          call. = FALSE )
  }

  hits <- lapply( X = site.list,
                  FUN = function(x){
                    findOverlaps( query = features,
                                  subject = x,
                                  minoverlap = min.overlap,
                                  type = "any",
                                  ignore.strand = TRUE ) } )

  frac.overlap <- lapply( X = hits,
                          FUN = function(x){
                            total.sites <- length( countSubjectHits(x) )
                            # overlapping <- length( countSubjectHits(x)[ countSubjectHits(x) != 0 ] )
                            overlapping <- length( unique( subjectHits(x) ) )
                            frac.overlap <- overlapping/total.sites
                            multioverlap <- length( countSubjectHits(x)[ countSubjectHits(x) > 1 ] )

                            if ( multioverlap != 0 ){
                              ol.check <- length( countSubjectHits(x)[ countSubjectHits(x) != 0 ] ) == overlapping

                              if( !isTRUE( ol.check ) ){
                                stop( "Error in counting feature overlaps.",
                                      call. = FALSE )
                                }
                              }

                            data.frame( total.sites = total.sites,
                                        overlapping.sites = overlapping,
                                        fraction.overlap = frac.overlap )
                            }
                          )

  frac.overlap <- do.call( rbind, frac.overlap )

  frac.overlap <- data.frame( sample = sample.names,
                              frac.overlap,
                              condition = factor( conditions, levels = condition.levels ) )

  feature.counts <- lapply( X = hits,
                            FUN = function(x){
                              counts <- c( countQueryHits(x) )
                              ids <- c( mcols( features )[ , names( mcols( features ) ) %in% id.col ] )
                              ids <- make.unique( ids, sep = ".dup" )
                              return( data.frame( ids, counts ) )
                              }
                            )

  feature.counts <- Map( cbind, feature.counts, sample = sample.names )

  # TO-DO: Rewrite this to remove tidyverse functions.
  count.mat <- as.matrix(
    do.call( rbind, c( feature.counts, make.row.names = FALSE ) ) |>
      tidyr::pivot_wider( names_from = sample, values_from = counts ) |>
      tibble::column_to_rownames( var = "ids" )
    )

  xint.obj <- SummarizedExperiment( assays = list( counts = count.mat ),
                                    colData = frac.overlap,
                                    rowRanges = features )

  xint.obj <- new( "xIntObject", xint.obj )

  return( xint.obj )
}
