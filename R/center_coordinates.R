center_coordinates <- function( site.list, current = 3, tsd = 5 ){

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

                  return( centered )
                  }
                )

  return( ll )

}
