#' Build Genome Mask
#'
#' Masks the target genome using a hidden Markov model according to genome fragmentation and sequencing methods. Assumes counts along genomic bins follow either a Poisson or negative binomial mixture model. Initial model parameters are estimated using one of fit_nb_mix() or fit_pois_mix(), assign_states(), and estimate_trans_mat(). Returns a GRanges object containing reduced genomic coordinates and their assigned state.
#'
#'@param aligned.fragments A GRanges object containing the coordinates of the genomic fragments extracted from the alignment BAM file.
#'@param genome.obj The BSgenome object for the target genome.
#'@param coarse.win.size The initial tile size with which to bin the genome. Defaults to 100.
#'@param reltol The relative tolerance for convergence. Used in the Baum-Welch algorithm to fit HMM parameters. Corresponds to the absolute relative difference between iterations. Defaults to 1E-10.
#'@param model The model to use in the EM algorithm. The options are <code>pois</code> and <code>nbinom</code>
#'@param ignore.chromosomes Character vector of chromosome names to ignore. These are not included in the HMM, and should therefore be removed from downstream analyses that depend on the HMM.
#'@param n.cores The number of cores to use during coarse mask trimming.
#'
#'@importFrom GenomeInfoDb seqlengths
#'@import GenomicRanges
#'@import HiddenMarkov
#'@import parallel
#'
#'@export
#'
mask_genome <- function( aligned.fragments,
                         genome.obj,
                         coarse.win.size = 100,
                         reltol = 1E-12,
                         model = c( "pois", "nbinom" ),
                         ignore.chromosomes = NULL,
                         n.cores = 1 ){

  n.states <- 2

  if( is( aligned.fragments, "list" ) || is( aligned.fragments, "GRangesList" ) ){
    frags <- aligned.fragments
  } else{
    if( is( aligned.fragments, "GRanges" ) ){
      frags <- list( aligned.fragments )
    } else{
      stop( "aligned.fragments must be a list of GRanges objects, a GRangesList, or a GRanges object.",
            call. = FALSE )
    }
  }

  cat( "Filtering chromosomes...", "\n\n" )

  if( !is.null( ignore.chromosomes ) ){
    frags <- mclapply( X = 1:length( frags ),
                       FUN = function(x) {
                         tmp <- frags[[x]]
                         og.len <- length( tmp )
                         tmp <- tmp[ !seqnames( tmp ) %in% ignore.chromosomes ]
                         new.len <- length( tmp )
                         cat( "Object", paste0( x, "/", length( frags ), ":" ),
                              "Filtered", og.len - new.len, "fragments out of", og.len, "original fragments",
                              paste0( "(",
                                      round( ( ( og.len - new.len ) / og.len ) * 100, 3 ), "%)."
                              ), "\n" )
                         return( tmp )
                         },
                       mc.cores = 12,
                       mc.preschedule = FALSE
                     )
    cat( "\n" )
  }

  nfrags <- lapply( X = frags,
                    FUN=function(x){
                      as.numeric( length(x) )
                      }
                    )

  nfrags <- sum( do.call( c, n.frags ) )

  cat( "Sorting fragments...", "\n\n" )

  frags <- mclapply( X = 1:length(frags),
                     FUN = function(x){
                       tmp <- frags[[x]]
                       if( isTRUE( is.unsorted( tmp ) ) ){
                         cat( "Sorting object",
                              paste0( x, " of ", length( frags ) ),
                              "\n" )
                         sort( tmp, ignore.strand = TRUE )
                       } else{
                         cat( "Object", x, "is sorted.",
                              "\n")
                       }
                     },
                     mc.cores = n.cores,
                     mc.preschedule = FALSE
                     )

  cat( "\n" )

  cat( "Tiling genome...", "\n\n" )

  sn <- lapply( X = frags,
                FUN = function(x){
                  unique( seqnames( x ) )
                  }
                )

  lvls <- unique( as.character( do.call( c, sn ) ) )

  sl <- seqlengths( genome.obj )[ lvls ]

  if( nfrags / sum( sl ) < 1 ){
    stop( "Insufficient number of fragments. Number of fragments / (actual) genome length should be >= 1." )
  }

  tiles <- tileGenome( seqlengths = sl, tilewidth = coarse.win.size, cut.last.tile.in.chrom = TRUE )
  mcols( tiles )$name <- paste0( "tile_", 1:length( tiles ) )
  ratio <- coarse.win.size / width( tiles )

  cat( "Counting tile overlaps...", "\n\n" )

  counts <- mclapply( X = 1:length(frags),
                      FUN = function(x){
                        cat( "Counting object", x, "overlaps...",
                             "\n" )
                        tmp <- frags[[x]]
                        countOverlaps( query = tiles,
                                       subject = tmp,
                                       minoverlap = 1,
                                       type = "any",
                                       ignore.strand = TRUE )
                        },
                      mc.cores = n.cores,
                      mc.preschedule = FALSE
                    )

  cat( "\n" )

  counts <- rowSums( do.call( cbind, counts ) )
  corr.counts <- round( counts * ratio )
  mcols( tiles )$count <- counts
  mcols( tiles )$corr.count <- corr.counts

  model <- match.arg( model )

  cat( "Deriving coarse mask...",
       "\n" )

  grls <- sapply( X = seqlevels( tiles ),
                  FUN = function( x ){

                    t1 <- Sys.time()

                    cat( "\n" )
                    cat( "Modeling chromosome ", gsub("chr", "", x), "\n" )
                    cat( "-----------------------","\n\n" )

                    tmp <- tiles[ seqnames( tiles ) == x ]
                    counts <- mcols( tmp )$corr.count

                    cat( length( tmp ), "tiles containing", sum( counts ), "counts.", "\n\n" )

                    cat( "Finding initial model parameters...", "\n\n" )

                    mmat <- switch( model,
                                    "pois" = fit_pois_mix( x = counts, size = 1E5 ),
                                    "nbinom" = fit_nb_mix( x = counts, size = 1E5 ) )
                    states <- assign_states( x = counts, mix.mat = mmat )
                    tmat <- estimate_trans_mat( x = states )

                    cat( "Training HMM...", "\n\n" )

                    if ( model == "nbinom" ){
                      init.params <- list( mu = mmat[, 2], size = mmat[, 3] )
                    } else{
                      init.params <- list( lambda = mmat[, 2] )
                    }

                    delta <- c( 0.5, 0.5 )

                    mod <- dthmm(
                      x = counts, Pi = tmat, delta = delta, distn = model,
                      pm = init.params, discrete = TRUE
                    )

                    ctrl <- bwcontrol(
                      maxiter = 1000, tol = reltol, prt = TRUE, posdiff = FALSE,
                      converge = expression( abs( diff / oldLL ) < tol )
                    )

                    mod.em <- bw( mod, control = ctrl )

                    ll <- logLik( mod.em, fortran = TRUE )

                    cat( "Final log-likelihood:",
                         ll,
                         "\n\n" )

                    cat( "Decoding hidden states...", "\n\n" )

                    mod.states <- Viterbi( mod.em )

                    mcols( tmp )$state <- mod.states

                    t2 <- Sys.time()

                    hmm.time <- round( as.double( difftime( t2, t1, units = "secs" ) ), 3 )

                    cat( hmm.time, "seconds elapsed.", "\n\n" )

                    tmp <- reduce( tmp[ mcols( tmp )$state == n.states ], ignore.strand = TRUE, min.gapwidth = 1L )

                    tot.mask <- sum( width( tmp ) )

                    cat( "Percent of chromosome coarse masked:",
                         paste0( round( ( tot.mask / sl[ x ] ) * 100, 3 ), "%" ), "\n\n" )

                    cat( "Done!", "\n\n" )

                    return( tmp )
                  },
                  simplify = FALSE )

  gr <- unlist( as( grls, "GRangesList" ) )

  gr <- sort( gr, ignore.strand = TRUE )

  cat( "\n", "Refining coarse mask.", "Utilizing", n.cores, "cores...", "\n\n" )

  cm.frags <- mclapply( X = frags,
                        FUN = function(x){
                          hits <- findOverlaps(query = gr,
                                               subject = x,
                                               minoverlap = 1L,
                                               type="any")
                          hits <- subjectHits( hits )
                          x[ hits ]
                          },
                        mc.cores = n.cores
                      )

  cm.frags <- sort( unlist( as( cm.frags, "GRangesList" ) ), ignore.strand = TRUE )

  trim <- refine_coarse_mask( coarse.mask = gr, fragments = cm.frags, n.cores = n.cores )

  cat( "Percent of each chromosome masked:", "\n" )

  invisible(
    lapply( X = names( sl ),
            FUN = function(x){
              val <- round( ( sum( width( trim[ seqnames( trim ) == x ] ) ) /
                                sl[ x ] ) * 100, 3 )
              cat( x, "\t", paste0( val, "%" ), "\n" )
            }
    )
  )

  cat( "\n" )

  cat( "Percent of genome masked:",
       paste0( round( ( sum( width( trim ) ) / sum( sl ) ) * 100, 3 ), "%" ), "\n\n" )

  return( trim )

}
