#' Make Fragments
#'
#' Make random fragments from random integration site positions and genomic cleavage sites
#'
#'@param insert.sites The generated random integration sites
#'@param frag.sites The genomic cleavage sites
#'@param random Boolean. Whether or not frag.sites were generated randomly (TRUE) or by simulated restriction digestion (FALSE)
#'@param mean Only required when random is TRUE. The target mean fragment size. 500 bp default.
#'@param sd Only required when random is TRUE. The target fragment size standard deviation. 250 bp default.
#'@param genome.obj The name of the genome object of interest.
#'
#'@import GenomeInfoDb
#'@import GenomicRanges
#'
#'@export
#'
make_fragments <- function( insert.sites,
                            frag.sites = NULL,
                            random = TRUE,
                            mean = 500,
                            sd = 250,
                            genome.obj ){

  genome.seqlengths <- seqlengths( genome.obj )

  if( !isTRUE( random ) ){
    if( is.null( frag.sites ) ){
      stop( "frag.sites cannot be NULL for non-random fragmentation.", call.=FALSE )
    }

    matches <- lapply( X = frag.sites,
                       FUN = function(x){
                         precede( x = insert.sites,
                                  subject = x,
                                  select = "all",
                                  ignore.strand = FALSE ) } )

    coordinates <- lapply( X=seq_along (matches ),
                           FUN=function(x){
                             matrix(c( start( insert.sites[ queryHits( matches[[ x ]] ) ] ),
                                       start( frag.sites[[ x ]][ subjectHits( matches[[ x ]] ) ] ) ),
                                    ncol = 2) } )

    frag.ranges <- lapply( X = seq_along( coordinates ),
                           FUN=function(x) {
                             names <- seqnames( insert.sites[ queryHits( matches[[ x ]] ) ] )
                             strand <- strand( insert.sites[ queryHits( matches[[ x ]] ) ] )
                             GRanges( seqnames = names,
                                      ranges = IRanges( start = apply( X = coordinates[[ x ]], MARGIN = 1, FUN = min ),
                                                        end = apply( X = coordinates[[ x ]], MARGIN = 1, FUN = max ) ),
                                      strand = strand ) } )

    frag.ranges <- do.call(c, frag.ranges)

    unmatched <- insert.sites[ !insert.sites %over% frag.ranges, ]

    unmatched.pos <- unmatched[ strand( unmatched ) == "+" ]
    pos.seqnames <- as.character( seqnames( unmatched.pos ) )
    end( unmatched.pos ) <- genome.seqlengths[ c( pos.seqnames ) ]

    unmatched.neg <- unmatched[ strand( unmatched ) == "-" ]
    neg.seqnames <- as.character( seqnames( unmatched.neg ) )
    start( unmatched.neg ) <- 1

    unmatched.frags <- do.call( c, list( unmatched.pos, unmatched.neg ) )

    frag.ranges <- do.call( c, list( frag.ranges, unmatched.frags ) )
  }

  else {
    lnorm.loc <- log( mean^2 / sqrt( sd^2 + mean^2 ) )
    lnorm.shape <- sqrt( log( 1 + ( sd^2 / mean^2 ) ) )

    plus.sites <- insert.sites[ strand( insert.sites ) == "+" ]
    plus.widths <- rlnorm( n = length(plus.sites), meanlog = lnorm.loc, sdlog = lnorm.shape)
    plus.sites <- GRanges( seqnames = seqnames(plus.sites),
                           ranges = IRanges(start = pmax( start( plus.sites ), 1 ),
                                            end = start( plus.sites ) + plus.widths ),
                           strand = "+" )

    minus.sites <- insert.sites[ strand( insert.sites ) == "-" ]
    minus.widths <- rlnorm( n = length( minus.sites ), meanlog = lnorm.loc, sdlog = lnorm.shape )
    minus.sites <- GRanges( seqnames = seqnames( minus.sites ),
                            ranges = IRanges(start = pmax( start( minus.sites ) - minus.widths, 1),
                                             end = end( minus.sites ) ),
                            strand = "-" )

    frag.ranges <- do.call( c, list( plus.sites, minus.sites ) )

    outliers <- bound_check( fragments = frag.ranges,
                             genome.obj = genome.obj,
                             include.lower = FALSE )

    if( length( outliers ) != 0 ){
      extreme.frags <- frag.ranges[ outliers ]
      extreme.seqnames <- as.character( seqnames( extreme.frags ) )
      end( extreme.frags ) <- genome.seqlengths[ c( extreme.seqnames ) ]
      frag.ranges <- frag.ranges[ -outliers ]
      frag.ranges <- do.call( c, list( frag.ranges, extreme.frags ) )
    }
  }

  frag.ranges <- sortSeqlevels( frag.ranges )
  frag.ranges <- sort( frag.ranges, ignore.strand = TRUE )

  return( frag.ranges )

}
