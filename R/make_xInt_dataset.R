#' Make xInt Dataset
#'
#' Creates a ranged summarized experiment object compatible with other xInt functions. Stores summary information for each sample in column data.
#'
#'@param site.list The list of integration site datasets
#'@param features The feature-set of interest
#'@param min.overlap The minimum amount of overlap to count two regions as overlapping
#'@param id.column The ID column in the features GRanges object. Default "name".
#'
#'@import GenomicRanges
#'@import SummarizedExperiment
#'@importFrom dplyr arrange
#'@importFrom tidyr pivot_wider
#'@importFrom tibble column_to_rownames
#'
#'@export
#'
make_xInt_dataset <- function( site.list, features, min.overlap = 1, id.col = "name" ){

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
                            overlapping <- length( countSubjectHits(x)[ countSubjectHits(x) != 0 ] )
                            # non.overlapping <- length( countSubjectHits(x)[  countSubjectHits(x) == 0 ] )
                            frac.overlap <- overlapping/total.sites
                            multioverlap <- length( countSubjectHits(x)[ countSubjectHits(x) > 1 ] )
                            if ( multioverlap != 0 ){
                               warning(cat( "Sites overlapping multiple features are present in the data.",
                                            "\n",
                                            "Use 'collapse_features()' to combined overlapping feature ranges if assessment of total feature integration is desired.",
                                       call. = FALSE )
                               )
                            }

                            data.frame( total.sites = total.sites,
                                        overlapping.sites = overlapping,
                                        # non.overlapping.sites = non.overlapping,
                                        fraction.overlap = frac.overlap)
                         })

  frac.overlap <- do.call( rbind, frac.overlap )

  feature.counts <- lapply( X = hits,
                            FUN = function(x){
                              # non.overlapping <- length( countSubjectHits(x)[ countSubjectHits(x) == 0 ] )
                              # counts <- c( non.overlapping, countQueryHits(x) )
                              counts <- c( countQueryHits(x) )
                              # ids <- c( "non.overlapping", mcols(features)[ , names( mcols(features)) %in% id.column ] )
                              ids <- c( mcols(features)[ , names( mcols( features) ) %in% id.col ] )
                              ids <- ids[ order( match( ids, id.col ) ) ]
                              return( data.frame( ids, counts ) )
                              } )

  feature.counts <- Map( cbind, feature.counts, sample = names( feature.counts ) )

  # Potential to-do: re-write the following to eliminate tidyverse functions. Maybe reshape + merge?
  count.mat <- as.matrix(
    do.call( rbind, c( feature.counts, make.row.names = FALSE ) ) |>
      dplyr::arrange( sample ) |>
      tidyr::pivot_wider( names_from = sample, values_from = counts ) |>
      tibble::column_to_rownames(var = "ids")
    )

  xint.obj <- SummarizedExperiment( assays = list( counts = count.mat ),
                                    colData = frac.overlap,
                                    rowRanges = features )

  return( xint.obj )
}
