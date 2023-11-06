#' Poisson Mixture Model
#'
#' Fits a Poisson mixture model to a vector of counts. Parameters are estimated by MLE. Generates initial mean estimates and mixing probabilities by k-means clustering. Returns an n x 2 matrix, where n is the number of distributions, column 1 holds the mixing proportions and column 2 holds the distribution means. The distribution with the highest mean value is returned in returned in row 1. The distribution with the smallest mean value is returned in row n.
#'
#'@param x The observed counts
#'@param n.dists The number of distributions in the mixture model.
#'@param sites The number of counts to sample from x.
#'
#'@export
#'
fit_pois_mix <- function( x, n.dists = 2, size = 1E5 ) {
  # Given a vector of counts, fit the distribution to a Poisson mixture distribution.
  # Generates initial mean and mixing probability estimates via k-means clustering.
  # Sorts final fitted values by decreasing mean values.
  # i.e. the largest mean is defined as state 1, second largest is defined as state 2, etc.
  # The output is used to assign states in assign_states().

  if ( n.dists < 2 ) {
    stop( "n.dists must be >= 2." )
    }

  # Define the negative log likelihood function
  nll <- function( params, x ) {
    lambdas <- params[ seq( 1, n.dists ) ]
    mps <- params[ seq( n.dists + 1, 2 * n.dists ) ]

    # Normalize the mixture probabilities
    mps <- mps / sum( mps )

    # Check for zero probabilities and adjust if necessary
    mps[ mps == 0 ] <- 1e-6

    # Compute the probabilities for each component of the mixture
    probs <- sapply( X = 1:n.dists,
                     FUN = function(i){
                       mps[i] * dpois( x, lambda = lambdas[i] )
                     } )

    # Compute the negative log likelihood
    calc.nll <- -sum( log( rowSums( probs, na.rm = TRUE ) ) )

    return( calc.nll )
  }

  max.tries <- 10
  tries <- 1

  while( tries <= max.tries ){
    tries <- tries + 1

    if ( length( x ) > size ){
      x <- sample( x = x, size = size, replace = FALSE )
    }

    # Generate initial parameter estimates using k-means clustering
    k <- kmeans( x, centers = n.dists )
    lambdas.init <- sort( k$centers )
    denom <- pmax( 1, 10^floor( log10( lambdas.init[1] ) ) )
    lambdas.init <- c( lambdas.init[1] / denom, lambdas.init[-1] )
    mps.init <- k$size / sum( k$size )

    # Generate random starting point for ps
    # ps_init <- runif(n.dists)

    # Combine lambdas and mps to create the initial parameter vector
    init.params <- c( lambdas.init, mps.init )

    # Minimize the negative log likelihood using the optim() function
    result <- tryCatch(
      {
        optim(
          fn = nll,
          par = init.params,
          method = "L-BFGS-B",
          lower = c( rep( 0, n.dists ), rep( 0, n.dists ) ),
          upper = c( rep( Inf, n.dists ), rep( 1, n.dists ) ),
          x = x
          )
        },
      error = function( e ) {
        NULL
      }
    )

    # Check if the optimization was successful
    if ( !is.null( result ) && result$convergence == 0 && !is.null( result$par ) ) {
      break
    } else{
      cat( "Optimization failed. Retrying...", "\n\n" )
    }
  }

  if (tries > max.tries) {
    stop( cat( "Optimization failed after", max.tries, "attempts.", "\n\n" ) )
  }

  # Extract the optimized parameters
  optimized.params <- result$par

  # Ensure mixture probabilities sum to 1
  optimized.params[ ( n.dists + 1 ):( 2 * n.dists ) ] <-
    optimized.params[ ( n.dists + 1 ):( 2 * n.dists ) ] /
    sum( optimized.params[ ( n.dists + 1 ):( 2 * n.dists ) ] )

  best.lambda <- optimized.params[ seq( 1, n.dists ) ]
  best.mp <- optimized.params[ seq( n.dists + 1, 2 * n.dists ) ]

  ind <- order( -best.lambda )

  best.lambda <- best.lambda[ind]
  best.mp <- best.mp[ind]

  final <- matrix( c( best.mp, best.lambda ), nrow = n.dists )
  colnames( final ) <- c( "mp", "lambda" )

  return( final )
}
