#' Create Genome Mask
#'
#' Enumerates regions of the genome that cannot be uniquely mapped
#' given a particular experimental approach (genome fragmentation method, read length, etc.).
#'
#'@param mapped.positions A GRanges object containing the positions of mapped fragments.
#'@param genome.obj The BSgenome object of interest.
#'@param ignore.chromosomes Chromosomes to be ignored in masking.
#'These chromosomes are not included in the final output.
#'@param n.cores The number of cores to use for parallel computing.
#'
#'@return A GRanges object containing the masked regions of the genome.
#'
#'@importFrom stats dpois rpois quantile
#'@import methods
#'@import changepoint
#'@import parallel
#'@import GenomicRanges
#'@importFrom IRanges IRanges %over%
#'
#'@export
#'
mask_genome <- function( mapped.positions,
                         genome.obj,
                         ignore.chromosomes = NULL,
                         n.cores = 1 ){

  # TO-DO: replace mclapply() with e.g. parLapply()

  t1 <- Sys.time()

  penalty <- "Hannan-Quinn"

  if( is( mapped.positions, "GRanges" ) ){
    frags <- mapped.positions
    } else{
      stop( "mapped.positions must be a GRanges object.",
            call. = FALSE )
    }

  message( "Filtering chromosomes...", "\n" )

  if( !is.null( ignore.chromosomes ) ){
    og.len <- length( frags )
    frags <- frags[ !seqnames( frags ) %in% ignore.chromosomes ]
    new.len <- length( frags )
    message( "Filtered ", og.len - new.len, " fragments out of ", og.len, " original fragments",
             " (",
             round( ( ( og.len - new.len ) / og.len ) * 100, 3 ), "%).",
             "\n" )
  } else{
      new.len <- length( frags )
    }

  lvls <- as.character( unique( seqnames( frags ) ) )

  if( n.cores < length( lvls ) ){
    warning( "n.cores is less than the number of chromosomes.
             This may cause forking issues due to insufficient memory.",
             "\n" )
    }

  nfrags <- new.len

  message( "Sorting fragments...", "\n" )

  if( is.unsorted( frags ) ){
    frags <- sort( frags, ignore.strand = TRUE )
    }

  sl <- seqlengths( genome.obj )[ lvls ]

  theo.mean <- ( nfrags / sum( sl ) )
  if( theo.mean < 1 ){
    stop( "Insufficient number of fragments. Number of fragments / genome length should be >= 1." )
  }

  mseg <- rle( rpois( 1E6, theo.mean ) == 0 )
  mseg <- quantile( mseg$lengths[ mseg$values == TRUE ], 0.99 )

  message( "Coarse masking genome...", "\n" )

  coarse.win.size <- mseg * 10

  coarse.tiles <- mclapply(
    X = names(sl),
    FUN = function(x){
      message( "Tiling chromosome ", x, " in ", coarse.win.size, " bp windows...",
               "\n" )
      if( ( sl[x] / coarse.win.size ) > 2^31-1 ){
        t <- GRanges( seqnames = NULL,
                      ranges = NULL,
                      strand = NULL )
        warning( "Chromosome", names( sl[x] ), "is longer than the maximum allowed vector length.",
                 "\n",
                 "Returning a NULL object.",
                 "\n" )
        } else{
          t <- tileGenome( seqlengths = sl[x], tilewidth = coarse.win.size, cut.last.tile.in.chrom = TRUE )
        }
      return(t)
      },
    mc.cores = n.cores,
    mc.preschedule = FALSE
   )

  message( "Counting coarse tile overlaps...", "\n" )

  coarse.tiles <- mclapply(
    X = seq_along(coarse.tiles),
    FUN = function(x){
      t <- coarse.tiles[[x]]
      counts <- countOverlaps( query = t,
                               subject = frags,
                               minoverlap = 1,
                               type = "any",
                               ignore.strand = TRUE )
     mcols( t )$count <- counts

     return( t )
     },
   mc.cores = n.cores,
   mc.preschedule = FALSE
   )

  message( "Extracting coarse mask...", "\n" )

  coarse.mask <- lapply(
    X = coarse.tiles,
    FUN = function(x){
      reduce( x[ mcols(x)$count == 0 ] )
    }
  )

  # The list elements in coarse.mask have no seqlevels in common.
  # suppressWarnings() is used here to suppress the warning telling us that.
  coarse.mask <- suppressWarnings( do.call( c, coarse.mask ) )

  message( "\n" )

  message( "Refining coarse mask...", "\n" )

  tiles <- mclapply(
    X = names(sl),
    FUN = function(x){
      message( "Tiling chromosome ", x, " in 1 bp windows...",
               "\n" )
      if( ( sl[x] ) > 2^31-1 ){
        t <- GRanges( seqnames = NULL,
                      ranges = NULL,
                      strand = NULL )
        warning( "Chromosome", names( sl[x] ), "is longer than the maximum allowed vector length.",
                 "\n",
                 "Returning a NULL object.",
                 "\n" )
      } else{
        t <- tileGenome( seqlengths = sl[x], tilewidth = 1, cut.last.tile.in.chrom = TRUE )
      }
      return(t)
    },
    mc.cores = n.cores,
    mc.preschedule = FALSE
  )

  tiles <- lapply(
    X = tiles,
    FUN = function(x){
      x[ !x %over% coarse.mask ]
    }
  )

  message( "Counting bp overlaps...", "\n" )

  tiles <- mclapply(
    X = seq_along(tiles),
    FUN = function(x){
      t <- tiles[[x]]
      counts <- countOverlaps( query = t,
                               subject = frags,
                               minoverlap = 1,
                               type = "any",
                               ignore.strand = TRUE )
      mcols( t )$count <- counts

      return( t )
    },
    mc.cores = n.cores,
    mc.preschedule = FALSE
  )

  message( "Finding changepoints...", "\n" )

  masked <- mclapply( X = tiles,
                      FUN = function(x){

                        t3 <- Sys.time()

                        n <- unique( as.character( seqnames( x ) ) )

                        message( "Finding chromosome", n, "changepoints...",
                                 "\n")

                        coi <- mcols(x)$count

                        chpt <- cpt.meanvar( data = coi,
                                             penalty = penalty,
                                             method = "PELT",
                                             test.stat = "Poisson",
                                             class = FALSE,
                                             minseglen = mseg )

                        coi.mat <- matrix( nrow = length( chpt ), ncol = 2 )

                        for ( i in seq_along( chpt ) ) {
                          if ( i == 1 ) {
                            coi.mat[i, 1] <- 1
                            coi.mat[i, 2] <- chpt[i]
                          } else {
                            coi.mat[i, 1] <- chpt[i - 1] + 1
                            coi.mat[i, 2] <- chpt[i]
                          }
                        }

                        coi.frac <- apply( X = coi.mat,
                                           MARGIN = 1,
                                           FUN = function(z) {
                                             v <- coi[ z[1]:z[2] ]
                                             f <- length( v[ v != 0 ] ) / length(v)
                                             f2 <- round( dpois( 0, theo.mean ) * length(v) ) / length(v)
                                             f < f2
                                             }
                                           )

                        regs <- coi.mat[ which( coi.frac ), ]

                        if( !is.matrix( regs ) ){
                          regs <- matrix( regs, ncol = 2 )
                        }
                        if( length( regs ) == 0 ){
                          gr <- GRanges( seqnames = NULL,
                                         ranges = NULL,
                                         strand = NULL )
                        } else{
                          starts <- regs[, 1]
                          stops <- regs[, 2]
                          a <- start( x[ starts ] )
                          b <- end( x[ stops ] )
                          gr <- GRanges( seqnames = n,
                                         ranges = IRanges( start = a, end = b ),
                                         strand = "*" )
                        }

                        t4 <- Sys.time()
                        xtime <- round( as.double( difftime( t4, t3, units = "secs" ) ), 3 )
                        message( "Done with chromosome ",
                                 n, "!",
                                 " ",
                                 xtime,
                                 " seconds elapsed...",
                                 "\n" )
                        return( gr )
                      },
                      mc.cores = n.cores,
                      mc.preschedule = FALSE )

  message( "\n" )
  message( "Combining coarse and fine masks...", "\n" )
  # Each element of 'masked' represents a different chromosome.
  # These elements, therefore, have no seqlevels in common.
  # suppressWarnings() is used here to suppress the warning telling us that.
  masked <- suppressWarnings( unlist( as( masked, "GRangesList" ) ) )
  masked <- sort( x = c( masked, coarse.mask ), ignore.strand = TRUE )
  masked <- reduce( x = masked, ignore.strand = TRUE )

  t2 <- Sys.time()
  tottime <- round( as.double( difftime( t2, t1, units = "secs" ) ), 3 )
  message( "Done! ", tottime, " total seconds elapsed.", "\n" )

  message( "Percent of each chromosome masked:", "\n" )

  invisible(
    lapply( X = names( sl ),
            FUN = function(x){
              val <- round( ( sum( width( masked[ seqnames( masked ) == x ] ) ) /
                                sl[ x ] ) * 100, 3 )
              message( x, "\t", val, "%", "\n" )
            }
    )
  )

  message( "\n" )

  message( "Percent of genome masked:",
           round( ( sum( width( masked ) ) / sum( sl ) ) * 100, 3 ), "%", "\n" )

  return( masked )

}
