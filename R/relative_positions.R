#' Relative Positions Along Features
#'
#' Calculates the relative position of each integration site along the features-of-interest. Returns either the relative positions for each site in each dataset or the number of sites within defined bins along the features.
#'
#'@param site.list The list or GRangesList containing the integration site coordinates.
#'@param features The features-of-interest. Must be a GRanges object.
#'@param bins The number of bins to generate across the features. Defaults to NULL, which returns the relative position of every site..
#'@param metadata The names of metadata columns to keep. Only used when the data are not binned.
#'@param ignore.strand Boolean. Whether or not to return relative positions relative to feature directionality. Defaults to FALSE (directionality is considered). This has no influence on overlap quantification.
#'
#'@return A list of data frames.
#'
#'@import GenomicRanges
#'@import S4Vectors
#'
#'@export
#'
relative_positions <- function( site.list,
                                features,
                                bins = NULL,
                                metadata = NULL ){

  check_sites( site.list )

  pos.ll <- lapply( X = site.list,
                    FUN = function(x){

                      if( all( width( x ) == 1 ) ){

                        gr <- x

                      } else{

                        warning( "The width of the provided integration sites was > 1.
                                 Defining the position of each integration site as the 5' end of the given coordinates.",
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

  if( !is.null( bins ) ){
    b <- seq(0, 1, length.out = bins + 1 )

    pos.bin <- lapply( X = pos.ll,
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
    return( pos.bin )
    } else{
      return( pos.ll )
    }
}
