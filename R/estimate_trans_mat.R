#' Estimate Transition Matrix
#'
#' Generate the initial estimate of the transition matrix for a hidden Markov model. The estimated vector of states can be generated with assign_states().
#'
#'@param x The estimated vector of states.
#'
#'
#'@export
#'
estimate_trans_mat <- function( x ) {
  # Estimates transition probabilities and creates a transition matrix for HMMs
  # Uses the states assigned in assign_states().
  # Generates reasonable starting values for the Baum-Welch algorithm.

  n.states <- 2

  # Enumerate the transitions present in the vector of states
  transpairs <- cbind( x[ -length( x ) ], x[ -1 ] )
  unique.states <- sort( unique( c( x, x + 1 ) ) )

  transprobs <- matrix( 0, nrow = n.states, ncol = n.states )

  # Count the number of times that state i transitions to state j
  # and create the transition matrix
  for ( i in 1:n.states ) {
    state.i <- unique.states[ i ]
    state.i.count <- length( which( x == state.i ) )

    for (j in 1:n.states) {
      state.j <- unique.states[ j ]

      state.j.count <- suppressWarnings(
        length( which(
          rowSums( transpairs ) == ( state.i + state.j ) & x == state.i ) )
        )

      if ( state.i.count == 0 ) {
        transprobs[ i,j ] <- 0
        } else {
          transprobs[ i,j ] <- state.j.count / state.i.count
        }
      }
    }

  colnames( transprobs ) <- paste( "State", 1:n.states )
  rownames( transprobs ) <- paste( "State", 1:n.states )

  return( transprobs )

}
