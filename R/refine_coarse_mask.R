refine_coarse_mask <- function( coarse.mask, fragments, n.cores = 1 ){

  t1 <- Sys.time()

  seqs <- as.character( unique( seqnames( coarse.mask ) ) )

  mtargs <- mclapply( X = seqs,
                      FUN = function(s){

                        cat( "Initiating chromosome", paste0( s, "..." ), "\n" )

                        f <- fragments[ seqnames( fragments ) == s ]
                        wt <- coarse.mask[ seqnames( coarse.mask ) == s ]
                        len <- length( wt )

                        tw <- lapply( X = 1:len,
                                      FUN = function(x){
                                        gr <- tile( x = coarse.mask[x],
                                                    width = 10 )
                                        gr <- unlist( gr )
                                        mcols( gr )$region <- x
                                        return( gr )
                                        }
                                      )

                        tw <- do.call( c, tw )

                        tw.counts <- countOverlaps( query = tw,
                                                    subject = f,
                                                    minoverlap = 1L,
                                                    type = "any",
                                                    ignore.strand = TRUE )

                        mcols( tw )$count <- tw.counts

                        tw <- lapply( X = unique( mcols( tw )$region ),
                                      FUN = function(x){

                                        region <<- x
                                        cat( "\r\033[K",
                                             paste0( "Chromosome ", s, ": "),
                                             "Trimming masked region", region,
                                             paste0( "of ", len, "..." ) )

                                        toi <- tw[ mcols( tw )$region == x ]
                                        coi <- mcols( toi )$count

                                        cpt <- cpt.np( data = coi,
                                                       minseglen = 1 )@cpts

                                        cpt <- cpt[ cpt != length( toi ) ]

                                        if( length( cpt ) == 0 ){
                                          GRanges( seqnames = NULL,
                                                   ranges = NULL,
                                                   strand = NULL )
                                          } else{
                                          cuts <- c( which( diff( cpt ) > 100 ), length( cpt ) )
                                          inits <- c( 1, cuts + 1 )
                                          inits <- inits[ inits <= length( cpt ) ]

                                          a <- start( toi[ cpt[ inits ] ] )
                                          b <- end( toi[ cpt[ cuts ] ] )

                                          GRanges( seqnames = unique( seqnames( toi ) ),
                                                   ranges = IRanges( start = a, end = b ),
                                                   strand = "*" )
                                          }
                                        }
                                      )

                        tw <- do.call( c, tw )

                        cat( "\n",
                             "Done with chromosome", paste0( s, "!" ),
                             "\n")

                        return( tw )

                        },
                        mc.cores = n.cores )

  mtargs <- do.call( c, mtargs )

  cat( "\n", "Subtracting targeted regions from coarse mask...", "\n\n" )

  new.mask <- subtract_mask( reference = coarse.mask, mask = mtargs )

  t2 <- Sys.time()

  tottime <- round( as.double( difftime( t2, t1, units = "secs" ) ), 3 )

  cat( "Done!", tottime, "seconds elapsed.", "\n\n" )

  return( new.mask )
}





