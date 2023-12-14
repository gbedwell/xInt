relative_positions <- function( site.list, features, names = NULL, metadata = NULL ){

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
                                        start = start( ranges( ol.feats ) ),
                                        end = end( ranges( ol.feats ) ) )

                      df$rel.position <- ( df$site - df$start ) / ( df$end - df$start )

                      if( !is.null( metadata ) ){
                        md <- mcols( ol.feats )[ names( mcols( features ) ) == metadata ]
                        df <- cbind( df, md )
                      }

                      return( df )
                      }
                    )

  if( is.null( names ) ){
    names( pos.ll ) <- names( site.list )
  } else{
    if( length( names ) != length( pos.ll ) ){
      stop( "Length of 'names' != length of 'site.list'.",
            "\n",
            "Verify that the provided names are correct." )
      } else{
        names( pos.ll ) <- names
      }
    }

  return( pos.ll )
}
