#' Global Integration Targeting Comparisons
#'
#' Compare global integration targeting preferences into a given feature-set between infection conditions. Comparisons are made using a two-sample t-test.
#'
#'@param xint.obj The xInt object of interest.
#'@param p.adjust.method The p-value adjustment method. Defaults to "BH". See <code>stats::p.adjust()</code> for more information.
#'
#'@return A data frame containing results of each comparison made.
#'
#'@importFrom stats t.test p.adjust
#'@import methods
#'
#'@export
#'
global_comparisons <- function( xint.obj, p.adjust.method = "BH" ){

  if( !validObject( xint.obj ) ){
    stop( "xint.obj is not a valid xIntObject.",
          call. = FALSE )
  }

  dat <- colData( xint.obj )

  lcpm <- log2( ( dat$overlapping.sites + 0.5 ) / ( dat$total.sites + 1 ) * 1E6 )
  names( lcpm ) <- dat$condition

  conditions <- levels( dat$condition )
  comparisons <- list()

  for (i in 1:( length( conditions ) - 1 ) ) {
    for ( j in ( i + 1 ):length( conditions ) ) {
      if ( conditions[i] != conditions[j] ){
        s1 <- lcpm[ names( lcpm ) == conditions[i] ]
        s2 <- lcpm[ names( lcpm ) == conditions[j] ]
        contrast.name <- paste( conditions[i], "-", conditions[j], sep = "" )
        tt <- t.test( x = s1,
                      y = s2,
                      alternative = "two.sided",
                      var.equal = TRUE )

        out <- data.frame(
          t = tt$statistic,
          df = tt$parameter,
          mu1 = tt$estimate[1],
          mu2 = tt$estimate[2],
          CI95.lower = tt$conf.int[1],
          CI95.upper = tt$conf.int[2],
          p.value = tt$p.value
        )

        comparisons[[contrast.name]] <- out
      }
    }
  }

  comparisons <- do.call( rbind, comparisons )
  comparisons$p.adj <- p.adjust( p = comparisons$p.value, method = p.adjust.method )

  return( comparisons )
}
