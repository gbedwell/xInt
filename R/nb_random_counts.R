#' Fit Random Counts to NB Distribution
#'
#' Holder text.
#'
#'@param sites A GRanges object containing the coordinates of the random genomic fragments.
#'@param mask A GRanges object containing the masked regions of the genome (e.g. the output of mask_genome()).
#'@param offset The nucleotide position of interest relative to the first nucleotide of the host sequence duplication. Defaults to 0, i.e. the first nucleotide of the duplication.
#'@param genome.obj The BSgenome object for the target genome.
#'@param size The number of observations to use in the fit. Defaults to 1E6.
#'@param win.size The tile size with which to bin the genome. Defaults to 100.
#'
#'@importFrom GenomeInfoDb seqlengths
#'@import GenomicRanges
#'
#'@export
#'
nb_random_counts <- function( sites,
                              mask,
                              offset = 0,
                              genome.obj,
                              size = 1E6,
                              win.size = 100 ){
  
  cat( "\n" )
  
  cat( "Checking range lengths...", "\n\n" )
  
  if( !all( width( sites ) == 1 ) ){
    
    cat( "Creating 1 bp site positions...", "\n\n" )
    
    plus <- sites[ strand( sites ) == "+" ]
    start( plus ) <- start( plus ) + offset
    end( plus ) <- start( plus )
    
    minus <- sites[ strand( sites ) == "-" ]
    start( minus ) <- end( minus ) - offset
    end( minus ) <- start( minus )
    
    sites <- sort( c( plus, minus ), ignore.strand = TRUE )
  }
  
  cat( "Tiling genome...", "\n\n" )
  
  lvls <- unique( seqnames( sites ) )
  
  sl <- seqlengths( genome.obj )[ lvls ]
  tiles <- tileGenome( seqlengths = sl, tilewidth = win.size, cut.last.tile.in.chrom = TRUE )
  mcols( tiles )$name <- paste0( "tile_", 1:length( tiles ) )
  
  cat( "Subtracting mask from tiles...", "\n\n" )
  
  tiles <- unlist( subtract( tiles, mask ) )
  
  ratio <- win.size / width( tiles )
  
  cat( "Counting tile overlaps...", "\n\n" )
  
  counts <- countOverlaps( query = tiles,
                           subject = sites,
                           minoverlap = 1,
                           type = "any",
                           ignore.strand = TRUE )
  
  tot.counts <- sum( counts )
  
  counts <- round( counts * ratio )
  
  if( length( counts ) < size ){
    size <- length( counts )
  }
  
  x <- sample( x = counts, size = size, replace = FALSE )
  
  # Define the negative log likelihood function
  nll <- function( params, x ) {
    mu <- params[1]
    phi <- params[2]
    
    # Compute the negative log likelihood
    calc.nll <- -sum( dnbinom( x, size = phi, mu = mu, log = TRUE ), na.rm = TRUE )
    
    return(calc.nll)
  }
  
  mu.init <- mean( counts )
  phi.init <- mu.init^2 / ( var( counts ) - mu.init )
  
  init.params <- c( mu.init, phi.init )
  
  cat( "Fitting NB distribution to observed counts...", "\n\n" )
  
  result <- optim( fn = nll,
                   par = init.params,
                   method = "L-BFGS-B",
                   lower = c( 0, 0 ),
                   upper = c( Inf, Inf ),
                   x = counts )
  
  df <- data.frame( win.size = win.size,
                    total = tot.counts,
                    mu = result$par[1],
                    phi = result$par[2] )
  
  return( df )
}



















