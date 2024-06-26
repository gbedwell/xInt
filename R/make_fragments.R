#' Make Fragments
#'
#' Make random fragments from random integration site positions and genomic cleavage sites
#'
#'@param insert.sites The generated random integration sites
#'@param frag.sites The positions in the genome where fragmentation occurs. For random fragmentation, this should be left NULL. If using restriction enzymes with pattern specificity, this should be the output of digest().
#'@param random Boolean. Whether or not frag.sites are to be generated randomly (TRUE) or by simulated restriction digestion (FALSE)
#'@param mean Only required when random is TRUE. The target mean fragment size. 500 bp default.
#'@param sd Only required when random is TRUE. The target fragment size standard deviation. 150 bp default.
#'@param genome.obj The BSgenome object of interest.
#'@param to.chr.ends Boolean. Whether or not to treat chromosome ends as capable of generating potentially mappable fragments. Defaults to TRUE.
#'
#'@return A GRanges object of fragment ranges.
#'
#'@importFrom stats rnbinom
#'@importFrom IRanges IRanges %over%
#'@import GenomeInfoDb
#'@import GenomicRanges
#'@import S4Vectors
#'
make_fragments <- function( insert.sites,
                            frag.sites = NULL,
                            random = TRUE,
                            mean = 500,
                            sd = 150,
                            genome.obj,
                            to.chr.ends = TRUE ){

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

    frag.ranges <- lapply( X = seq_along( matches ),
                           FUN = function(x){

                             is <- insert.sites[ queryHits( matches[[x]] ) ]
                             fs <- frag.sites[[x]][ subjectHits( matches[[x]] ) ]

                             plus.ind <- which( as.character( strand(is) ) == "+" )
                             minus.ind <- which( as.character( strand(is) ) == "-" )

                             plus <- GRanges( seqnames = seqnames( is[ plus.ind ]),
                                              ranges = IRanges( start = start( is[ plus.ind ] ),
                                                                end = end( fs[ plus.ind ] ) ),
                                              strand = strand( is[ plus.ind ] )
                             )

                             minus <- GRanges( seqnames = seqnames( is[ minus.ind ] ),
                                               ranges = IRanges( start = start( fs[ minus.ind ] ),
                                                                 end = end( is[ minus.ind ] ) ),
                                               strand = strand( is[ minus.ind ] )
                             )

                             gr <- c( plus, minus )
                             gr <- sort( gr, ignore.strand = TRUE )

                             return( gr )
                           }
    )

    frag.ranges <- do.call(c, frag.ranges)

    if( to.chr.ends ){

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
  } else {
    # Old: use log-normal distribution to sample random reads.
    # Superseded by use of negative binomial (below).
    # lnorm.loc <- log( mean^2 / sqrt( sd^2 + mean^2 ) )
    # lnorm.shape <- sqrt( log( 1 + ( sd^2 / mean^2 ) ) )

    disp <- mean^2 / ( sd^2 - mean )

    plus.sites <- insert.sites[ strand( insert.sites ) == "+" ]
    # plus.widths <- rlnorm( n = length(plus.sites), meanlog = lnorm.loc, sdlog = lnorm.shape)
    plus.widths <- rnbinom( n = length( plus.sites ), size = disp, mu = mean )
    plus.sites <- GRanges( seqnames = seqnames(plus.sites),
                           ranges = IRanges(start = start( plus.sites ),
                                            end = start( plus.sites ) + plus.widths ),
                           strand = "+" )

    minus.sites <- insert.sites[ strand( insert.sites ) == "-" ]
    # minus.widths <- rlnorm( n = length( minus.sites ), meanlog = lnorm.loc, sdlog = lnorm.shape )
    minus.widths <- rnbinom( n = length( minus.sites ), size = disp, mu = mean )
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
      end( extreme.frags ) <- genome.seqlengths[ extreme.seqnames ]
      frag.ranges <- frag.ranges[ -outliers ]
      frag.ranges <- do.call( c, list( frag.ranges, extreme.frags ) )
    }
  }

  frag.ranges <- sortSeqlevels( frag.ranges )
  frag.ranges <- sort( frag.ranges, ignore.strand = TRUE )

  return( frag.ranges )

}
