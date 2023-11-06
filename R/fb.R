#' Forward-Backward Algorithm
#'
#' A modified version of the forwardback() function in the HiddenMarkov package. This modified function was first written by Reed A. Cartwright (https://gist.github.com/reedacartwright/7b5339862c9abd388e10b101bfed6f07). The modifications help prevent numeric underflow by centering probability calculations around probmax. See the HiddenMarkov manual for more information on the stock function (https://cran.r-project.org/web/packages/HiddenMarkov/HiddenMarkov.pdf).
#'
#'@param x is a vector of length n containing the univariate observed process. Alternatively, x could be specified as NULL, meaning that the data will be added later (e.g. simulated).
#'@param Pi is the m × m transition probability matrix of the homogeneous hidden Markov chain
#'@param delta is the marginal probability distribution of the m hidden states at the first time point
#'@param distn is a character string with the abbreviated distribution name. Distributions provided by the package are Beta ("beta"), Binomial ("binom"), Exponential ("exp"), GammaDist ("gamma"), Lognormal ("lnorm"), Logistic ("logis"), Normal ("norm"), and Poisson ("pois"). See topic Mstep, Section “Modifications and Extensions”, to extend to other distributions.
#'@param pm is a list object containing the (Markov dependent) parameter values associated with the distribution of the observed process.
#'@param pn is a list object containing the observation dependent parameter values associated with the distribution of the observed process.
#'@param fortran logical, if TRUE (default) use the Fortran code, else use the R code.
#'
#'@export
#'
fb <- function (x, Pi, delta, distn, pm, pn = NULL, fortran = TRUE) {
  # Modified from https://gist.github.com/reedacartwright/7b5339862c9abd388e10b101bfed6f07
  m <- nrow(Pi)
  n <- length(x)
  dfunc <- HiddenMarkov:::makedensity(distn)
  prob <- matrix(as.double(0), nrow=n, ncol=m)
  for (k in 1:m)
    prob[,k] <- dfunc(x=x, HiddenMarkov:::getj(pm, k), pn, log=TRUE)
  probmax = apply(prob,1,max)
  prob = exp(prob-probmax)

  #   forward probabilities alpha_ij
  phi <- as.double(delta)
  logalpha <- matrix(as.double(rep(0, m*n)), nrow=n)
  lscale <- as.double(0)
  if (fortran!=TRUE){
    #  loop1 using R code
    for (i in 1:n){
      if (i > 1) phi <- phi %*% Pi
      phi <- phi*prob[i,]
      sumphi <- sum(phi)
      phi <- phi/sumphi
      lscale <- lscale + log(sumphi)
      logalpha[i,] <- log(phi) + lscale
    }
    LL <- lscale
  } else{
    if (!is.double(Pi)) stop("Pi is not double precision")
    if (!is.double(prob)) stop("prob is not double precision")
    memory0 <- rep(as.double(0), m)
    loop1 <- .Fortran("loop1", m, n, phi, prob, Pi, logalpha,
                      lscale, memory0, PACKAGE="HiddenMarkov")
    logalpha <- loop1[[6]]
    LL <- loop1[[7]]
  }

  #   backward probabilities beta_ij
  logbeta <- matrix(as.double(rep(0, m*n)), nrow=n)
  phi <- as.double(rep(1/m, m))
  lscale <- as.double(log(m))
  if (fortran!=TRUE){
    #  loop2 using R code
    for (i in seq(n-1, 1, -1)){
      phi <- Pi %*% (prob[i+1,]*phi)
      logbeta[i,] <- log(phi) + lscale
      sumphi <- sum(phi)
      phi <- phi/sumphi
      lscale <- lscale + log(sumphi)
    }
  } else{
    memory0 <- rep(as.double(0), m)
    loop2 <- .Fortran("loop2", m, n, phi, prob, Pi, logbeta,
                      lscale, memory0, PACKAGE="HiddenMarkov")
    logbeta <- loop2[[6]]
  }

  logalpha = logalpha + cumsum(probmax)
  logbeta = logbeta + rev(cumsum(rev(c(probmax[-1],0))))
  LL = LL + sum(probmax)

  return(list(logalpha=logalpha, logbeta=logbeta, LL=LL))
}
