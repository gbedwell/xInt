center_coordinates <- function( site.list, current = 3, tsd = 5, genome.obj ){

  ll <- lapply( X = site.list,
                FUN = function(x){
                  center <- ( tsd + 1 ) / 2
                  diff <- center - current
                  a <- floor( diff )
                  b <- ceiling( diff )

                  plus <- x[ strand( x ) == "+" ]
                  end( plus ) <- start( plus ) + b
                  start( plus ) <- start( plus ) + a

                  minus <- x[ strand( x ) == "-" ]
                  start( minus ) <- end( minus ) - b
                  end( minus ) <- end( minus ) - a

                  centered <- sort( c( plus, minus ), ignore.strand = TRUE )

                  sn <- as.character( seqnames( centered ) )
                  sl <- seqlengths( genome.obj )

                  if( any( end( centered ) > ( sl[ sn ] - ceiling( center ) ) ) ||
                      any( start( centered ) < floor( center ) ) ){

                    warning( "1 or more centered coordinates are out of bounds. ",
                             "Returning the problematic ranges." )

                    problems <- c( centered[ which( end( centered ) > ( sl[ sn ] - ceiling( center ) ) ) ],
                                   centered[ which( start( centered ) < floor( center ) ) ] )

                    problems <- sort( problems, ignore.strand = TRUE )

                    } else{

                    return( centered )

                    }
                  }
                )

  return( ll )

}
