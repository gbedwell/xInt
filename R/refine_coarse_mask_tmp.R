#' Refine Coarse Mask
#'
#' Refines a coarse genomic mask by identifying non-masked sub-regions within the coarsely masked regions via change point estimation. This is generally not intended to be run alone, but is called within mask_genome().
#'
#'@param coarse.mask The GRanges object containing the coarsely masked genomic regions.
#'@param fragments The GRanges object containing the simulated random fragments/sites.
#'@param n.cores The number of cores to use in parallelization.
#'
#'@export
#'
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
                                        gr <- tile( x = wt[x],
                                                    width = 1 )
                                        gr <- unlist( gr )
                                        mcols( gr )$region <- x
                                        return( gr )
                                      }
                        )

                        tw <- unname( do.call( c, tw ) )

                        tw.counts <- countOverlaps( query = tw,
                                                    subject = f,
                                                    minoverlap = 1L,
                                                    type = "any",
                                                    ignore.strand = TRUE )

                        mcols( tw )$count <- tw.counts

                        tw <- lapply( X = unique( mcols( tw )$region ),
                                      FUN = function(x) {

                                        region <<- x

                                        cat( "\r\033[K",
                                             paste0( "Chromosome ", s, ": " ),
                                             "Trimming masked region", region,
                                             paste0( "of ", len, "..." ) )

                                        toi <- unname( tw[ mcols( tw )$region == x ] )
                                        coi <- unname( mcols( toi )$count )

                                        chpt <- cpt.meanvar( data = coi,
                                                             penalty = "MBIC",
                                                             method = "PELT",
                                                             test.stat = "Poisson",
                                                             class = FALSE,
                                                             minseglen = 10 )

                                        if ( length( chpt ) == 1 && chpt == length( coi ) ) {

                                          tmp.frac <- length( coi[ coi != 0 ] ) / length( coi )

                                          if( tmp.frac >= 0.05 ){

                                            a <- start( toi[ 1 ] )
                                            b <- end( toi[ length( toi ) ] )

                                            GRanges( seqnames = unique( seqnames( toi ) ),
                                                     ranges = IRanges( start = a, end = b ),
                                                     strand = "*" )

                                          } else{

                                            GRanges( seqnames = NULL,
                                                     ranges = NULL,
                                                     strand = NULL )

                                          }
                                        } else {

                                          coi.mat <- matrix( nrow = length( chpt ), ncol = 2 )

                                          for ( i in 1:length( chpt ) ) {
                                            if ( i == 1 ) {

                                              coi.mat[i, 1] <- 1
                                              coi.mat[i, 2] <- chpt[i]

                                            } else {

                                              coi.mat[i, 1] <- chpt[i - 1] + 1
                                              coi.mat[i, 2] <- chpt[i]

                                            }
                                          }

                                          coi.mean <- apply( X = coi.mat,
                                                             MARGIN = 1,
                                                             FUN = function(x) {
                                                               mean( coi[x[1]:x[2]] )
                                                             }
                                          )

                                          coi.frac <- apply( X = coi.mat,
                                                             MARGIN = 1,
                                                             FUN = function(x) {
                                                               v <- coi[ x[1]:x[2] ]
                                                               length( v[ v != 0 ] ) / length(v)
                                                             }
                                          )

                                          regs <- coi.mat[ which( coi.frac >= 0.01 ), ]

                                          if( !is.matrix( regs ) ){
                                            regs <- matrix( regs, ncol = 2 )
                                          }

                                          if( length( regs ) == 0 ){

                                            GRanges( seqnames = NULL,
                                                     ranges = NULL,
                                                     strand = NULL )

                                          } else{

                                            starts <- regs[, 1]
                                            stops <- regs[, 2]

                                            a <- start( toi[ starts ] )
                                            b <- end( toi[ stops ] )

                                            GRanges( seqnames = unique( seqnames( toi ) ),
                                                     ranges = IRanges( start = a, end = b ),
                                                     strand = "*" )
                                          }
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

  mtargs <- reduce( x = mtargs, ignore.strand = TRUE )

  cat( "\n", "Subtracting targeted regions from coarse mask...", "\n\n" )

  new.mask <- subtract_mask( reference = coarse.mask, mask = mtargs )

  t2 <- Sys.time()

  tottime <- round( as.double( difftime( t2, t1, units = "secs" ) ), 3 )

  cat( "Done!", tottime, "seconds elapsed.", "\n\n" )

  return( new.mask )
}





