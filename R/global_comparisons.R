#' Global Integration Targeting Comparisons
#'
#' Compare global integration targeting preferences into a given feature-set between infection conditions.
#' Comparisons are made using a two-sample t-test.
#'
#'@param xint.obj The xInt object of interest.
#'@param p.adjust.method The p-value adjustment method.
#'Defaults to "BH". See <code>stats::p.adjust()</code> for more information.
#'
#'@return A data frame containing the results of each comparison.
#'
#'@examples
#'data(xobj)
#'global_comparisons(xint.obj = xobj)
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

  # Implement same log-cpm calculation in limma.
  # Ensures that log2(0) is never taken
  # and that overlapping.sites/total.sites is always < 1.
  lcpm <- log2( ( dat$overlapping.sites + 0.5 ) / ( dat$total.sites + 1 ) * 1E6 )
  names( lcpm ) <- dat$condition

  conditions <- levels( dat$condition )

  cond.lens <- vapply(
    X = conditions,
    FUN = function(x){
      length( dat$condition[ dat$condition == x ] )
      },
    FUN.VALUE = numeric(1)
    )

  multi.rep <- names( cond.lens[ cond.lens > 1 ] )
  single.rep <- names( cond.lens[ cond.lens == 1 ] )

  if( length( multi.rep ) == 0 ){
    stop( "No replicate datasets are provided.",
          call. = FALSE )
  }

  comparisons <- list()

  for( i in seq( 1, ( length( multi.rep ) - 1 ) ) ){
    for( j in seq( ( i + 1 ), length( multi.rep ) ) ){
      if( multi.rep[i] != multi.rep[j] ){
        s1 <- lcpm[ names( lcpm ) == multi.rep[i] ]
        s2 <- lcpm[ names( lcpm ) == multi.rep[j] ]
        contrast.name <- paste( multi.rep[i], "-", multi.rep[j], sep = "" )
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

  if( length( single.rep ) != 0 ){
    for( i in seq( 1, length( multi.rep ) ) ){
      for( j in seq( 1, length( single.rep ) ) ){
        if( multi.rep[i] != single.rep[j] ){
          s1 <- lcpm[ names( lcpm ) == multi.rep[i] ]
          s2 <- lcpm[ names( lcpm ) == single.rep[j] ]
          contrast.name <- paste( multi.rep[i], "-", single.rep[j], sep = "" )
          tt <- t.test( x = s1,
                        mu = s2,
                        alternative = "two.sided",
                        var.equal = TRUE )

          out <- data.frame(
            t = tt$statistic,
            df = tt$parameter,
            mu1 = tt$estimate[1],
            mu2 = tt$null.value,
            CI95.lower = tt$conf.int[1],
            CI95.upper = tt$conf.int[2],
            p.value = tt$p.value
          )

          comparisons[[contrast.name]] <- out
        }
      }
    }
    warning( "The following datasets cannot be compared: ",
             paste( single.rep, collapse = ", " ), ".",
             call. = FALSE )
  }

  comparisons <- do.call( rbind, comparisons )
  comparisons$p.adj <- p.adjust( p = comparisons$p.value, method = p.adjust.method )

  return( comparisons )
}
