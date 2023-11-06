#' Expectation Algorithm
#'
#' A modified version of the Estep() function in the HiddenMarkov package. The function uses a modified forward-backward algorithm (fb(), instead of the HiddenMarkov forwardback()). Parameter definitions are taken directly from the HiddenMarkov manual (https://cran.r-project.org/web/packages/HiddenMarkov/HiddenMarkov.pdf).
#'
#'@param x is a vector of length n containing the univariate observed process. Alternatively, x could be specified as NULL, meaning that the data will be added later (e.g. simulated).
#'@param Pi is the m × m transition probability matrix of the homogeneous hidden Markov chain
#'@param delta is the marginal probability distribution of the m hidden states at the first time point
#'@param distn is a character string with the abbreviated distribution name. Distributions provided by the package are Beta ("beta"), Binomial ("binom"), Exponential ("exp"), GammaDist ("gamma"), Lognormal ("lnorm"), Logistic ("logis"), Normal ("norm"), and Poisson ("pois"). See topic Mstep, Section “Modifications and Extensions”, to extend to other distributions.
#'@param pm is a list object containing the (Markov dependent) parameter values associated with the distribution of the observed process.
#'
#'@export
#'
estep <- function (x, Pi, delta, distn, pm, pn = NULL){
    dfunc <- HiddenMarkov:::makedensity(distn)
    m <- nrow(Pi)
    n <- length(x)
    y <- fb(x, Pi, delta, distn, pm, pn)
    logbeta <- y$logbeta
    logalpha <- y$logalpha
    LL <- y$LL
    u <- exp(logalpha + logbeta - LL)
    v <- array(NA, dim = c(n - 1, m, m))
    for (k in 1:m) {
      logprob <- dfunc(x=x[-1], HiddenMarkov:::getj(pm, k),
                       HiddenMarkov:::getj(pn, -1), log=TRUE)
      logPi <- matrix(log(Pi[, k]), byrow = TRUE, nrow = n -
                        1, ncol = m)
      logPbeta <- matrix(logprob + logbeta[-1, k],
                         byrow = FALSE, nrow = n - 1, ncol = m)
      v[, , k] <- logPi + logalpha[-n, ] + logPbeta - LL
    }
    v <- exp(v)
    return(list(u = u, v = v, LL = LL))
  }
