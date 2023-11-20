#' Classic Proportion Comparisons
#'
#' Performs global pairwise comparisons on an xInt object using classic tests for equal proportions.
#'
#'@param xint.obj The xInt object containing the relevant data.
#'@param comparison The extent of the pairwise comparisons. <code>all</code> will compare every dataset to all other datasets. <code>pooled</code> pools the information by condition and performs pairwise comparisons across conditions.
#'@param p.adj The method by which to correct p-values for multiple comparisons. Can be NULL. Defaults to "BH". See p.adjust() documentation for more details.
#'
#'@export
#'
classic_prop_tests <- function( xint.obj, comparison = c( "all", "pooled" ), p.adj = "BH" ){

  # Potential to-do: implement tests to return log(p) in case of extremely small p-values.

  comparison = match.arg( comparison )

  g.test <- function( x ){
    ex <- outer( rowSums( x ), colSums( x ) )/sum( x )
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
    labels <- unique( unfactor( colData( xint.obj )$condition ) )
    total <- with( colData( xint.obj ), tapply( total.sites, condition, sum ) )
    n.in <- with( colData( xint.obj ), tapply( overlapping.sites, condition, sum ) )
    n.out <- total - n.in
    frac <- n.in / total
    combos <- combn( labels, 2 )
  }

  pairwise <- sapply(X = 1:ncol( combos ),
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
                     })

  rownames( pairwise ) <- c( "Fisher", "Chi.Sq", "G", "OR", "h" )

  or <- pairwise[ 4, ]
  pairh <- pairwise[ 5, ]

  if (!is.null( p.adj ) && ncol( pairwise ) > 1 ){
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

  if( any( pairwise[, 1:3 ] == 0 ) ){
    warning( "Some comparisons returned p = 0.",
             "\n",
             paste0( "This is because the associated probability is < ", signif( .Machine$double.xmin, 4 ), ", the smallest non-zero normalized floating-point number." ),
             "\n",
             paste0( "For practical purposes, reporting these p-values as e.g., < ", signif( .Machine$double.eps, 2 ), " is recommended." ),
             "\n" )
  }

  return( pairwise )
}
