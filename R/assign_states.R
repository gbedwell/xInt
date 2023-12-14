#' Assign Initial States
#'
#' Initial HMM state assignments for a vector of counts based on defined thresholds. States are assigned such that the state above the highest threshold is State 1, the state above the next highest threshold is State 2, etc.
#'
#'@param x The observed count vector.
#'@param mix.mat The parameter matrix output by fit_nb_mix() or fit_pois_mix().
#'@param thresholds A vector of thresholds that define the different states. An alternative to supplying mix.mat.
#'
#'@export
#'
assign_states <- function( x, mix.mat = NULL, thresholds = NULL ) {
  # Assign states to the vector of counts to be input into the HMM.
  # Written to be compatible with fit_pois_mix() and fit_nb_mix() output.
  # Thresholds can also be defined manually.
  # Number of thresholds should be the n.states - 1.
  # This is used to approximate transition probabilities in estimate_trans_mat().

  n.states <- 2

  if ( !is.null( mix.mat ) ){
    n.states <- nrow( mix.mat )
    thresholds <- mix.mat[ 2:nrow( mix.mat ), 2 ]
  } else{
    n.states <- 2
    if ( is.null( thresholds ) ){
      stop( "One of mix.mat or thresholds must be defined.", call. = FALSE )
    }
    if ( length( thresholds ) != n.states - 1 ){
      stop( "The number of defined thresholds must be n.states-1.", call. = FALSE )
    }
  }

  thresholds <- sort( thresholds, decreasing = TRUE )

  if ( any( thresholds %% 1 != 0 ) ){
    thresholds <- floor( thresholds )
  }

  states <- rep( 1, length( x ) )

  for ( i in 2:n.states ) {
    states[ x < thresholds[ i - 1 ] ] <- i
  }

  return( states )

}
