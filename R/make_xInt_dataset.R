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
#'@importFrom tidyr pivot_wider
#'@importFrom tibble column_to_rownames
#'
#'@export
#'
make_xInt_dataset <- function( site.list,
                               features,
                               sample.names = NULL,
                               conditions,
                               condition.levels,
                               min.overlap = 1,
                               id.col = "name" ){

  if( is.null( sample.names ) ){
    sample.names <- names( site.list )
    if( is.null( sample.names ) ){
      stop( "Sample names must be provided or the elements in site.list must be named.",
            call. = FALSE )
    }
  }

  if( length( sample.names ) != length( site.list ) ){
    stop( "The number of sample names does not match the number of integration site datasets.",
          call. = FALSE )
  }

  if( length( sample.names ) != length( conditions ) ){
    stop( "The number of sample names does not match the number of conditions.",
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

                               warning( "Sites overlapping multiple features are present in the data.",
                                        "\n",
                                        "The summarized sample (column) data will not be affected.",
                                        "\n",
                                        call. = FALSE )

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
                              ids <- c( mcols( features )[ , names( mcols( features) ) %in% id.col ] )
                              ids <- make.unique( ids, sep = ".dup" )
                              return( data.frame( ids, counts ) )
                              }
                            )

  feature.counts <- Map( cbind, feature.counts, sample = sample.names )

  count.mat <- as.matrix(
    do.call( rbind, c( feature.counts, make.row.names = FALSE ) ) |>
      tidyr::pivot_wider( names_from = sample, values_from = counts ) |>
      tibble::column_to_rownames(var = "ids")
    )

  xint.obj <- SummarizedExperiment( assays = list( counts = count.mat ),
                                    colData = frac.overlap,
                                    rowRanges = features )

  if( all( colnames( assay( xint.obj ) ) == rownames( colData( xint.obj ) ) ) ){
    return( xint.obj )
    } else{
    stop( "Sample naming issue.",
          call. = FALSE )
  }
}
