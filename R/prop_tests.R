#' Proportion Comparisons
#'
#' Performs global pairwise comparisons on an xInt object using classic tests for equal proportions.
#'
#'@param xint.obj The xIntObject containing the relevant data.
#'@param comparison The extent of the pairwise comparisons.
#'<code>all</code> will compare every dataset to all other datasets.
#'This can take a long time.
#'<code>pooled</code> pools the information by condition and performs pairwise comparisons across conditions.
#'@param p.adj The method by which to correct p-values for multiple comparisons.
#'Defaults to "BH". See p.adjust() documentation for more details.
#'
#'@examples
#'data(xobj)
#'prop_tests(xint.obj = xobj, comparison="pooled")
#'
#'@importFrom utils combn
#'@importFrom stats chisq.test fisher.test pchisq p.adjust
#'@import methods
#'@import SummarizedExperiment
#'
#'@return A matrix of p-values and effect sizes for each pairwise comparison performed.
#'
#'@export
#'
prop_tests <- function( xint.obj, comparison = c( "all", "pooled" ), p.adj = "BH" ){

  if( !validObject( xint.obj ) ){
    stop( "xint.obj is not a valid xIntObject.",
          call. = FALSE )
  }

  comparison = match.arg( comparison )

  g.test <- function( x ){
    ex <- outer( rowSums( x ), colSums( x ) ) / sum( x )
    G <- 2 * sum( x * log( x / ex ) )
    pchisq( q = G, df = 1, lower.tail = FALSE, log = FALSE)
  }

  if( comparison == "all" ){
    labels <- row.names( colData( xint.obj ) )
    frac <- xint.obj$fraction.overlap
    total <- xint.obj$total.sites
    n.in <- xint.obj$overlapping.sites
    n.out <- xint.obj$total.sites - xint.obj$overlapping.sites
    combos <- combn( labels, 2 )
    names( frac ) <- names( total ) <- names( n.in ) <- names( n.out ) <- labels
  } else{
    labels <- unique( as.character( colData( xint.obj )$condition ) )
    total <- with( colData( xint.obj ), tapply( total.sites, condition, sum ) )
    n.in <- with( colData( xint.obj ), tapply( overlapping.sites, condition, sum ) )
    n.out <- total - n.in
    frac <- n.in / total
    combos <- combn( labels, 2 )
  }

  pairwise <- vapply(X = seq_len( ncol( combos ) ),
                     FUN = function(x){
                       col <- combos[,x]
                       mat <- matrix( c( n.in[col[1]],
                                         n.in[col[2]],
                                         n.out[col[1]],
                                         n.out[col[2]] ),
                                      nrow = 2)
                       ft <- fisher.test(x = mat, alternative = "two.sided")$p.value
                       chi <- chisq.test( x = mat )$p.value
                       g <- g.test( x = mat )
                       or <- ( n.in[col[1]] / n.out[col[1]] ) / ( n.in[col[2]] / n.out[col[2]] )
                       p1 <- 2 * asin( sqrt( frac[ col[1] ] ) )
                       p2 <- 2 * asin( sqrt( frac[ col[2] ] ) )
                       pairh <- abs( p1 - p2 )
                       rbind( ft, chi, g, or, pairh )
                       },
                     FUN.VALUE = numeric( 5 )
                     )

  rownames( pairwise ) <- c( "Fisher", "Chi.Sq", "G", "OR", "h" )

  or <- pairwise[ 4, ]
  pairh <- pairwise[ 5, ]

  if ( ncol( pairwise ) > 1 ){
    pairwise <- apply(X = pairwise[1:3,],
                      MARGIN = 1,
                      FUN = function(x){
                        p.adjust(p = x, method = p.adj)
                      })
    pairwise <- cbind( pairwise, or, pairh )
  } else{
    pairwise <- t( pairwise )
  }

  colnames( pairwise )[4:5] <- c( "OR", "h" )

  rownames( pairwise ) <- apply(X = combos, MARGIN = 2, FUN = function(x){ paste0( x[1], "-", x[2] ) })

  # This should handle NA and NaN values.
  # OR the functions above should not return NaN...
  if( any( pairwise[, 1:3 ] == 0 ) ){
    warning( "Some comparisons returned p = 0.",
             "\n",
             "This is because the associated probability is < ",
             signif( .Machine$double.xmin, 4 ),
             ", the smallest non-zero normalized floating-point number.",
             "\n",
             "For practical purposes, reporting these p-values as e.g., < ",
             signif( .Machine$double.eps, 2 ),
             " is recommended.",
             "\n",
             call. = FALSE )
  }

  return( pairwise )
}
