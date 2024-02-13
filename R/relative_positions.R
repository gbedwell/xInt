#' Meta-feature Analysis
#'
#' Calculates the relative positions of each site along the features-of-interest. Can return the relative positions for each site in each dataset, or the number of sites within defined bins along the features.
#'
#'@param site.list The list or GRangesList containing the mapped site coordinates.
#'@param features The features-of-interest. Must be a GRanges object.
#'@param metadata The names of metadata columns to keep. Only works when <code>bins=NULL</code>. That is, when the relative positions of every site is returned.
#'@param bins The number of bins to generate across the features. Defaults to NULL (i.e., the relative position of every site is returned).
#'
#'@export
#'
relative_positions <- function( site.list,
                                features,
                                metadata = NULL,
                                bins = NULL ){

  pos.ll <- lapply( X = site.list,
                    FUN = function(z){

                      if( all( width( z ) == 1 ) ){

                        gr <- z

                      } else{

                        plus <- z[ strand( z ) == "+" ]
                        end( plus ) <- start( plus )

                        minus <- z[ strand( z ) == "-" ]
                        start( minus ) <- end( minus )

                        gr <- sort( c( plus, minus ), ignore.strand = TRUE )

                        }

                      ol <- findOverlaps( query = gr,
                                          subject = features,
                                          minoverlap = 1L,
                                          type = "any" )

                      ol.sites <- gr[ queryHits( ol ) ]
                      ol.feats <- features[ subjectHits( ol ) ]


                      df <- data.frame( site = start( ranges( ol.sites ) ),
                                        feature.start = start( ranges( ol.feats ) ),
                                        feature.end = end( ranges( ol.feats ) ) )

                      df$rel.position <- ( df$site - df$feature.start ) / ( df$feature.end - df$feature.start )

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
