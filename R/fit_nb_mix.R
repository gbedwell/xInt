#' Negative Binomial Mixture Model
#'
#' Fits a negative binomial mixture model to a vector of counts. Parameters are estimated by MLE. Generates initial mean estimates and mixing probabilities by k-means clustering. Returns an n x 2 matrix, where n is the number of distributions, column 1 holds the mixing proportions and column 2 holds the distribution means. The distribution with the highest mean value is returned in returned in row 1. The distribution with the smallest mean value is returned in row n.
#'
#'@param x The observed counts
#'@param size The number of counts to sample from x.
#'
#'@export
#'
fit_nb_mix <- function(x, size = 1E5) {
  # Given a vector of counts, fit the distribution to a negative binomial mixture distribution.
  # Generates initial mean and mixing probability estimates via k-means clustering.
  # Sorts final fitted values by decreasing mean values.
  # i.e. the largest mean is defined as state 1, second largest is defined as state 2, etc.
  # The output is used to assign states in assign_states().

  n.dists <- 2

  # Define the negative log likelihood function
  nll <- function( params, x ) {
    mus <- params[ seq( 1, n.dists ) ]
    phis <- params[ seq( n.dists + 1, 2 * n.dists ) ]
    mps <- params[ seq (2 * n.dists + 1, 3 * n.dists ) ]

    # Normalize the mixture probabilities
    mps <- mps / sum( mps )

    # Check for zero probabilities and adjust if necessary
    mps[ mps == 0 ] <- 1e-6

    # Compute the probabilities for each component of the mixture
    probs <- sapply(X = 1:n.dists, FUN = function(i) {
      mps[i] * dnbinom( x, size = phis[i], mu = mus[i] )
    })

    # Compute the negative log likelihood
    calc.nll <- -sum( log( rowSums( probs, na.rm = TRUE ) ) )

    return( calc.nll )
  }

  max.tries <- 5
  tries <- 1

  while ( tries <= max.tries ) {
    tries <- tries + 1

    if ( length( x ) > size ) {
      samp <- sample( x = x, size = size, replace = FALSE )
    }

    dens <- stats::density( x = samp,
                            from = 0,
                            to = max( samp ),
                            n = length( 0:max( samp ) ) )$y

    inf.pt <- which( diff( sign( diff( dens ) ) ) == 2 ) + 1

    p1 <- samp[ which( samp >= inf.pt[1] ) ]
    p2 <- samp[ which( samp < inf.pt[1] ) ]

    mus.init <- c( mean( p1 ), mean( p2 ) )
    vars.init <- c( var( p1 ), var( p2 ) )
    phis.init <- mus.init^2 / ( vars.init - mus.init )
    phis.init[ phis.init <= 0 ] <- 1
    mps.init <- c( ( length( p1 ) / length( samp ) ), ( length( p2 ) / length( samp ) ) )

    # Combine mus, phis, and mps to create the initial parameter vector
    init.params <- c( mus.init, phis.init, mps.init )


    # Minimize the negative log likelihood using the optim() function
    result <- tryCatch(
      {
        optim(
          fn = nll,
          par = init.params,
          method = "L-BFGS-B",
          lower = c(rep(0, n.dists), rep(0, n.dists), rep(0, n.dists)),
          upper = c(rep(Inf, n.dists), rep(Inf, n.dists), rep(1, n.dists)),
          x = samp
        )
      },
      error = function(e) {
        NULL
      }
    )

    # Check if the optimization was successful
    if (!is.null(result) && result$convergence == 0 && !is.null(result$par) && !any(is.na(result$par))) {
      break
    } else {
      cat("Optimization failed. Retrying...", "\n\n")
    }
  }

  if (tries > max.tries) {
    warning( "Optimization failed after ", max.tries, " attempts.", "\n",
             "Using initial parameter estimates instead.", "\n\n")
    optimized.params <- init.params
  }

  # Extract the optimized parameters
  optimized.params <- result$par

  # Ensure mixture probabilities sum to 1
  optimized.params[(2 * n.dists + 1):(3 * n.dists)] <-
    optimized.params[(2 * n.dists + 1):(3 * n.dists)] /
    sum(optimized.params[(2 * n.dists + 1):(3 * n.dists)])

  best.mu <- optimized.params[seq(1, n.dists)]
  best.phi <- optimized.params[seq(n.dists + 1, 2 * n.dists)]
  best.mp <- optimized.params[seq(2 * n.dists + 1, 3 * n.dists)]

  ind <- order( -best.mu )

  best.mu <- best.mu[ind]
  best.phi <- best.phi[ind]
  best.mp <- best.mp[ind]

  final <- matrix( c( best.mp, best.mu, best.phi ), nrow = n.dists )
  colnames(final) <- c( "mp", "mu", "phi" )

  return(final)
}
