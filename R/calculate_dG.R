#' Calculate B-to-A deltaG
#'
#' Calculates the free energy change upon B to A transition of DNA.
#' Relies on trimeric dG values calculated by Tolstorukov, et al. (DOI: 10.1016/S0006-3495(01)75973-5).
#' Utilizes a 3 bp sliding window with a 1 bp step across the sequence of interest.
#'
#'@param site.list The list or GRangesList containing the mapped site coordinates.
#'@param step.len The target length of the returned vector of dG values.
#'@param genome.obj The genome object of interest.
#'@param current.start The position in the target site duplication currently described by start.
#'This is used for centering the site coordinates.
#'@param tsd The total length of the target site duplication.
#'This is used for centering the site coordinates.
#'@param collapse Boolean. Whether or not to collapse the output into a single data.frame. Defaults to TRUE.
#'
#'@return A list of data frames containing the sequence position and the corresponding dG value for each element in site.list.
#'
#'@examples
#'library(BSgenome.Hsapiens.UCSC.hs1)
#'data(sites2)
#'calculate_dG( site.list = sites2, genome.obj = Hsapiens )
#'
#'@import ggplot2
#'@import GenomicRanges
#'@import Biostrings
#'
#'@export
#'
calculate_dG <- function( site.list, step.len = 50, genome.obj,
                          current.start = 1, tsd = 5, collapse = TRUE ){

  if( !validObject( site.list ) ){
    stop( "site.list is not a valid SiteListObject.",
          call. = FALSE )
  }

  if( step.len %% 2 != 0 ){
    warning( "step.len cannot be odd. Subtracting 1 from step.len.",
             call. = FALSE )
    step.len <- step.len - 1
  }

  seq.len <- step.len + 2

  expanded.sites <- expand_coordinates( site.list = site.list,
                                        seq.len = seq.len,
                                        genome.obj = genome.obj,
                                        current.start = current.start,
                                        tsd = tsd )

  seqs <- getSeq( x = genome.obj,
                  names = as( expanded.sites, "GRangesList" ),
                  as.character = FALSE )

  trimer <-  c( "AAA", "AAC", "AAG", "AAT", "ATC", "ACA",
                "ACG", "AGA", "AGC", "AGG", "ATA", "CAA",
                "CAC", "CAG", "CCA", "CCG", "CGA", "CGC",
                "CTA", "GAA", "GAC", "GCA", "GCC", "GGA",
                "GTA", "TAA", "TCA", "ACC", "ACT", "ATG",
                "CCC", "CTC" )

  vals <- c( 1.17, 1.22, 1.57, 1.39, 1.11, 0.52,
             0.57, 0.95, 0.48, 0.64, 0.63, 0.52,
             0.67, 0.46, 0.44, 0.42, 0.56,0.64,
             0.66, 0.72, 0.49, 0.76, 0.54, 0.56,
             0.67, NA, 0.58, 0.3, 0.3, 0.22,
             0.26, 0.21 )

  vals = c( vals, vals )

  trimer.rev <- as.character( reverseComplement( DNAStringSet( trimer ) ) )

  names( vals ) <- c( trimer, trimer.rev )

  len <- width( expanded.sites[[1]][1] )

  if( len %% 2 == 1 ){
    len <- len - 1
  }

  dG <- lapply( X = seqs,
                FUN = function(x){
                  out <- calc_dG_internal(x, vals, len)

                  center <- ( tsd + 1 ) / 2

                  if( tsd %% 2 == 1 ){
                    mp <- length(out) / 2
                    zp <- mp - ( center - 1 )
                  } else{
                    mp <- length(out) / 2
                    zp <- mp - ( floor(center) - 1 )
                  }
                  dat <- data.frame( position = ( 0:( length(out) - 1 ) ) - ( zp ), dG = out )
                  if( tsd %% 2 == 1 ){
                    dat <- dat[ -c(1), ]
                  }
                  return( dat )
                  }
                )


  if( isTRUE( collapse ) ){
    dG <- Map( function(x, y) { cbind( x, dataset = y ) }, dG, names( dG ) )
    dG <- do.call( rbind, dG )
    rownames( dG ) <- NULL
  }

  return( dG )
}
