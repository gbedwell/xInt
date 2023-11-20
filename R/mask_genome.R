#' Build Genome Mask
#'
#' Masks the target genome using a hidden Markov model according to genome fragmentation and sequencing methods. Assumes counts along genomic bins follow either a Poisson or negative binomial mixture model. Initial model parameters are estimated using one of fit_nb_mix() or fit_pois_mix(), assign_states(), and estimate_trans_mat(). Returns a GRanges object containing reduced genomic coordinates and their assigned state.
#'
#'@param aligned.fragments A GRanges object containing the coordinates of the genomic fragments extracted from the alignment BAM file.
#'@param generated.fragments A GRanges object containing the coordinates of the genomic fragments generated via <code>make_fragments()</code>. These represent ground-truth fragments.
#'@param genome.obj The BSgenome object for the target genome.
#'@param coarse.win.size The initial tile size with which to bin the genome. Defaults to 100.
#'@param refine.win.size The window size at which to refine the coarse mask.
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
                         generated.fragments = NULL,
                         verify.frags = FALSE,
                         genome.obj,
                         coarse.win.size = 100,
                         refine.win.size = 10,
                         reltol = 1E-12,
                         model = c( "pois", "nbinom" ),
                         ignore.chromosomes = NULL,
                         n.cores = 1 ){

  n.states <- 2

  if( isTRUE( verify.frags ) && is.null( generated.fragments ) ){
    stop( "generated.fragments cannot be NULL when verify.frags = TRUE.",
          call. = FALSE )
  }

  cat( "\n" )

  cat( "Filtering chromosomes...", "\n\n" )

  og.len <- length( aligned.fragments )

  if ( !is.null( ignore.chromosomes ) ){
    aligned.fragments <- aligned.fragments[ !seqnames( aligned.fragments ) %in% ignore.chromosomes ]
  }

  new.len <- length( aligned.fragments )

  cat( "Filtered", og.len - new.len, "fragments out of", og.len, "original fragments",
       paste0( "(",
               round( ( ( og.len - new.len ) / og.len ) * 100, 3 ), "%)."
       ), "\n\n" )

  if( isTRUE( verify.frags ) ){
    cat( "Verifying fragments...", "\n\n" )

    frags <- verify_aligned_fragments( generated.fragments = generated.fragments,
                                       aligned.fragments = aligned.fragments )

    cat( "Keeping", length( frags ), "fragments out of",
         length( aligned.fragments ), "aligned fragments",
         paste0( "(",
                 round( ( length( frags ) / length( aligned.fragments ) ) * 100, 3 ), "%)."
         ),
         "\n\n" )

  } else{
    frags <- aligned.fragments
  }

  # cat( "Getting fragment positions...", "\n\n" )

  # f.len <- length( frags )
  # f.ind <- sort( sample( x = 1:f.len, size = ceiling( f.len / 2 ), replace = FALSE ) )
  #
  # five <- frags[ f.ind ]
  # fp <- five[ strand( five ) == "+" ]
  # start( fp ) <- start( fp )
  # end( fp ) <- start( fp )
  #
  # fm <- five[ strand( five ) == "-" ]
  # start( fm ) <- end( fm )
  # end( fm ) <- start( fm )
  #
  # three <- frags[ -f.ind ]
  # tp <- three[ strand( three ) == "+" ]
  # start( tp ) <- start( tp )
  # end( tp ) <- start( tp )
  #
  # tm <- three[ strand( three ) == "-" ]
  # start( tm ) <- end( tm )
  # end( tm ) <- start( tm )

  # offset <- 0
  #
  # plus <- frags[ strand( frags ) == "+" ]
  # start( plus ) <- start( plus )
  # end( plus ) <- start( plus )
  #
  # minus <- frags[ strand( frags ) == "-" ]
  # start( minus ) <- end( minus )
  # end( minus ) <- start( minus )

  if( !isTRUE( is.unsorted( frags ) ) ){
    cat( "Sorting fragments...", "\n\n" )
    frags <- sort( frags, ignore.strand = TRUE )
  }

  # frags <- sort( c( plus, minus ), ignore.strand = TRUE )

  # frags <- sort( c( fp, fm, tp, tm ), ignore.strand = TRUE )

  cat( "Tiling genome...", "\n\n" )

  lvls <- unique( seqnames( frags ) )

  sl <- seqlengths( genome.obj )[ lvls ]
  tiles <- tileGenome( seqlengths = sl, tilewidth = coarse.win.size, cut.last.tile.in.chrom = TRUE )
  mcols( tiles )$name <- paste0( "tile_", 1:length( tiles ) )
  ratio <- coarse.win.size / width( tiles )

  cat( "Counting tile overlaps...", "\n\n" )

  counts <- countOverlaps( query = tiles,
                           subject = frags,
                           minoverlap = 1,
                           type = "any",
                           ignore.strand = TRUE )

  corr.counts <- round( counts * ratio )
  mcols( tiles )$count <- counts
  mcols( tiles )$corr.count <- corr.counts

  model <- match.arg( model )

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
                                    "pois" = fit_pois_mix( x = counts, n.dists = n.states, size = 1E5 ),
                                    "nbinom" = fit_nb_mix( x = counts, n.dists = n.states, size = 1E5 ) )
                    states <- assign_states( x = counts, mix.mat = mmat )
                    tmat <- estimate_trans_mat( x = states, n.states = n.states )

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

                    k <- length( mod.em$pm ) * length( mod.em$pm[[ 1 ]] )

                    cat( "AIC:",
                         -2 * ll + 2 * k,
                         "\n\n" )

                    cat( "BIC:",
                         -2 * ll + k * log( length( tmp ) ),
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

  # gr <- reduce( gr, min.gapwidth = 1L )

  gr <- sort( gr, ignore.strand = TRUE )

  cat( "\n", "Refining coarse mask.", "Utilizing", n.cores, "cores...", "\n\n" )

  if( ceiling( sum( sl ) / length( frags ) ) > refine.win.size ){
    warning( "refine.win.size is too small given data depth.",
             "\n",
             paste0( "Using minimum allowable refine.win.size of ", ceiling( sum( sl ) / length( frags ) ), "." ),
             call. = FALSE )
    refine.win.size <- ceiling( sum( sl ) / length( frags ) )
  }

  trim <- refine_coarse_mask( coarse.mask = gr, fragments = frags, win.size = refine.win.size, n.cores = n.cores )

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
