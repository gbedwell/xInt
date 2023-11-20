#' Negative Binomial Maximization Function
#'
#' Maximization algorithm for negative binomial distributions. Modified from a function written by Reed A. Cartwright (https://gist.github.com/reedacartwright/7b5339862c9abd388e10b101bfed6f07).
#'
#'@param x is a vector of length n containing the univariate observed process. Alternatively, x could be specified as NULL, meaning that the data will be added later (e.g. simulated).
#'@param cond is an object created by estep.
#'@param pm is a list object containing the (Markov dependent) parameter values associated with the distribution of the observed process.
#'@param pn is a list object containing the observation dependent parameter values associated with the distribution of the observed process.
#'
#'@export
#'
Mstep.nbinom <- function(x, cond, pm, pn) {
  # Modified from https://gist.github.com/reedacartwright/7b5339862c9abd388e10b101bfed6f07
  w = cond$u
  mu = c()
  size = c()
  for(i in seq_len(ncol(w))) {
    ww = w[,i]
    A = weighted.mean(x, ww)
    r = pm$size[i]

    if( r <= 0 ){
      warning( "r <= 0, setting r = 1E-5", call. = FALSE )
      r <- 1E-5
    }

    # function for calculating the gradient
    g = function(r) {
      u = digamma(x+r)-digamma(r)
      u = u-log1p(A/r)
      as.vector(ww %*% u)
    }

    # function for calculating the hessian
    h = function(r) {
      v = trigamma(x+r)-trigamma(r)
      v = v + A/r/(A+r)
      as.vector(ww %*% v)
    }

    # Find the root of the gradient
    tryCatch( {
      o = uniroot(g, c(r*0.1, r*1.9), extendInt = "yes")
      size = c(size, o$root)
      mu = c(mu, A)
    }, error = function(e) {
      # Handle the error here (e.g., print a message or take alternative action)
      cat("An error occurred in uniroot: ", conditionMessage(e), "\n")
      # You can also expand the interval here if you want to retry with a wider range
      o = uniroot(g, c(r*1E-4, r*200), extendInt = "yes")
      size = c(size, o$root)
      mu = c(mu, A)
    })
  }

  # Construct return value
  return( list(mu=mu,size=size) )
}
