#' Refine Coarse Mask
#'
#' Refines a coarse genomic mask by identifying non-masked sub-regions within the coarsely masked regions via change point estimation. This is generally not intended to be run alone, but is called within mask_genome().
#'
#'@param coarse.mask The GRanges object containing the coarsely masked genomic regions.
#'@param fragments The GRanges object containing the simulated random fragments/sites.
#'@param win.size The window size used to bin the coarse masked regions. Defaults to 10.
#'@param n.cores The number of cores to use in parallelization.
#'
#'@export
#'
refine_coarse_mask <- function( coarse.mask, fragments, win.size = 10, n.cores = 1 ){

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
                                                    width = win.size )
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

                                        toi <- unname( tw[ mcols( tw )$region == x ] )
                                        coi <- unname( mcols( toi )$count )

                                        chpt <- cpt.meanvar( data = coi,
                                                             penalty = "MBIC",
                                                             method = "PELT",
                                                             test.stat = "Poisson",
                                                             class = FALSE,
                                                             minseglen = 1 )

                                        if( length( chpt ) == 1 && chpt == length( coi ) ){
                                          GRanges( seqnames = NULL,
                                                   ranges = NULL,
                                                   strand = NULL )
                                          } else{

                                            if( 1 %in% chpt ){

                                              coi.grp <- cut( 1:length( coi ), breaks = chpt )
                                              coi.splt <- split( coi, coi.grp )
                                              coi.splt <- c( `(1,1]` = coi[1], coi.splt )

                                            } else{

                                              coi.grp <- cut( 1:length( coi ), breaks = c( 1, chpt ) )
                                              coi.splt <- split( coi, coi.grp )
                                              coi.splt[[1]] <- c( coi[1], coi.splt[[1]] )

                                            }

                                            coi.mean <- unlist(
                                              lapply( X = coi.splt,
                                                      FUN = function(x){
                                                        mean(x)
                                                        }
                                                      )
                                              )

                                            regs <- lapply( X = names( coi.splt )[ which( coi.mean > 0.1 ) ],
                                                            FUN = function(x){
                                                              as.numeric(
                                                                scan(
                                                                  text = gsub( pattern = "\\(",
                                                                               replacement = "",
                                                                               x = gsub(
                                                                                 pattern = "\\]",
                                                                                 replacement = "",
                                                                                 x = x ) ),
                                                                  sep = ",",
                                                                  quiet = TRUE )
                                                                )
                                                              }
                                                            )

                                            regs <- do.call( rbind, regs )

                                            starts <- regs[,1]
                                            stops <- regs[,2]

                                            a <- start( toi[ starts ] )
                                            b <- end( toi[ stops ] )

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





