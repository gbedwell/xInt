#' Meta-feature analysis
#'
#' Calculates the relative position of each integration site along the features-of-interest.
#' Returns either the relative positions for each site in each dataset or the number of sites within defined
#' bins along the features.
#'
#'@param site.list The list or GRangesList containing IS coordinates.
#'@param features The features-of-interest. Must be a GRanges object.
#'@param bins The number of bins to generate across the features.
#'Defaults to NULL, which returns the relative position of every site.
#'<bins=0</code> behaves like NULL.
#'@param metadata The names of metadata columns to keep.
#'Only used when the data are not binned.
#'@param collapse Boolean. Whether or not to collapse the output into a single data frame.
#'Defaults to TRUE.
#'
#'@return A data frame or a list of data frames.
#'
#'@examples
#'data(sites2)
#'data(xobj)
#'feats <- rowRanges(xobj)
#'
#'metafeature(site.list = sites2,
#'            features = feats)
#'
#'metafeature(site.list = sites2,
#'            features = feats,
#'            bins = 10)
#'
#'@import GenomicRanges
#'@import S4Vectors
#'
#'@export
#'
metafeature <- function( site.list,
                         features,
                         bins = NULL,
                         metadata = NULL,
                         collapse = TRUE ){

  if( !validObject( site.list ) ){
    stop( "site.list is not a valid SiteListObject.",
          call. = FALSE )
  }

  pos.ll <- lapply( X = site.list,
                    FUN = function(x){

                      if( all( width( x ) == 1 ) ){

                        gr <- x

                      } else{

                        warning( "The width of the provided integration site coordinates is > 1.
                                 Defining the position of each integration site as the 5'
                                 end of the given coordinates.",
                                 call. = FALSE )

                        plus <- x[ strand( x ) == "+" ]
                        end( plus ) <- start( plus )

                        minus <- x[ strand( x ) == "-" ]
                        start( minus ) <- end( minus )

                        gr <- sort( c( plus, minus ), ignore.strand = TRUE )

                        }

                      ol <- findOverlaps( query = gr,
                                          subject = features,
                                          minoverlap = 1L,
                                          type = "any",
                                          ignore.strand = TRUE )

                      ol.sites <- gr[ queryHits( ol ) ]
                      ol.feats <- features[ subjectHits( ol ) ]


                      df <- data.frame( chr = as.character( seqnames( ol.sites ) ),
                                        site = start( ranges( ol.sites ) ),
                                        site.strand = as.character( strand( ol.sites ) ),
                                        feature.start = start( ranges( ol.feats ) ),
                                        feature.end = end( ranges( ol.feats ) ),
                                        feature.strand = as.character( strand( ol.feats ) ) )

                      df$rel.position <- ( df$site - df$feature.start ) / ( df$feature.end - df$feature.start )

                      df[ df$feature.strand == "-", ]$rel.position <- 1 - df[ df$feature.strand == "-", ]$rel.position

                      if( !is.null( metadata ) ){
                        md <- mcols( ol.feats )[ names( mcols( features ) ) == metadata ]
                        df <- cbind( df, md )
                      }

                      return( df )
                      }
                    )

  if( !is.null( bins ) & bins != 0 ){
    b <- seq(0, 1, length.out = bins + 1 )

    pos.ll <- lapply(
      X = pos.ll,
      FUN = function(x){
        v <- x$rel.position
        tt <- table(cut(v, breaks = b, include.lowest = TRUE))
        percentiles <- b
        percentiles <- percentiles[ percentiles != 0 ] * 100
        dat <- data.frame( percentiles = percentiles )
        dat <- cbind( dat, data.frame( tt ) )
        dat <- dat[,c(1,3)]
        colnames(dat) <- c( "percentile", "frequency" )
        dat$fraction <- dat$frequency / sum( dat$frequency )
        return(dat)
        }
      )
    }

  if( isTRUE( collapse ) ){
    pos.ll <- Map( function(x, y) { cbind( x, dataset = y ) }, pos.ll, names( pos.ll ) )
    pos.ll <- do.call( rbind, pos.ll )
    rownames( pos.ll ) <- NULL
  }

  return( pos.ll )
}
