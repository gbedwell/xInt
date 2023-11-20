relative_locations <- function( site.list, features, names = NULL, metadata = NULL, n.cores = 1 ){

  if( is.null( metadata ) ){
    metadata = names( mcols( features ) )
  }

  loc.ll <- mclapply( X = 1:length( site.list ),
                      FUN = function(z){
                        is <- site.list[[z]]

                        ol <- findOverlaps( query = is,
                                            subject = features,
                                            minoverlap = 1L,
                                            type = "any" )

                        sh <- unique( subjectHits( ol ) )

                        ll <- lapply( X = 1:length( sh ) ,
                                      FUN = function(x){
                                        tmp.ol <- ol[ subjectHits( ol ) == sh[ x ] ]
                                        foi <- features[ unique( subjectHits( tmp.ol ) ) ]
                                        sites <- data.frame( site = start( is[ queryHits( tmp.ol ) ] ) )
                                        df <- data.frame( site = sites,
                                                          ranges( foi ),
                                                          mcols( foi )[ metadata ] )
                                        df <- df[ !names( df ) %in% c( "width" ) ]
                                        df$rel.location <- (df$site - df$start ) / ( df$end - df$start )
                                        return( df )
                                        }
                                      )
                        ll <- do.call( rbind, ll )
                      },
                      mc.cores = n.cores
                      )

  names( loc.ll ) <- names( site.list )

  return( loc.ll )
}
