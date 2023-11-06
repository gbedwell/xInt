#' Baum-Welch algorithm
#'
#' A modified version of the BaumWelch.dthmm() function in the HiddenMarkov package. This function is modified to use the absolute relative difference as the convergence criteria instead of absolute difference. This is intended to make the function more capable of dealing with relatively flat likelihood functions. See the original package manual for more information on its inner-workings (https://cran.r-project.org/web/packages/HiddenMarkov/HiddenMarkov.pdf).
#'
#'@param object An object of class "dthmm", "mmglm0", "mmglm1", "mmglmlong1", or "mmpp".
#'@param control A list of control settings for the iterative process. These can be changed by using the function bwcontrol.
#'
#'@export
#'
bw <- function (object, control = bwcontrol(), ...){
  # modifies BaumWelch.dthmm() in the HiddenMarkov package
  # HiddenMarkov package manual: https://cran.r-project.org/web/packages/HiddenMarkov/HiddenMarkov.pdf
  # Original function source code: https://rdrr.io/cran/HiddenMarkov/src/R/BaumWelch.dthmm.R
  # Redefines the first iteration's oldLL from -Inf to -2^31-1 -- an arbitrarily large number and the largest storable integer in R.
  # This enables setting the convergence criterion to the relative difference instead of the absolute difference
  # i.e. abs( ( LL(n + 1) - LL(n) ) / LL(n) ) instead of abs( LL(n + 1) - LL(n) )
  # This is used in case the likelihood function is relatively flat (especially in the case of very large N)
  # In this case, the absolute difference between iterations might never reach small absolute difference values
  # Using relative difference enables consistent criterion use across all data and models
  # Also adds a relative difference output when printing iteration values.

  x <- object$x
  Pi <- object$Pi
  delta <- object$delta
  distn <- object$distn
  pm <- object$pm
  tol <- control$tol
  if (distn[1]!="glm"){
    Mstep <- parse(text=paste("Mstep.", distn,
                              "(x, cond, pm, object$pn)", sep=""))
  } else{
    Mstep <- parse(text=paste("Mstep.glm",
                              "(x, cond, pm, object$pn, distn[2], distn[3])", sep=""))
  }
  m <- nrow(Pi)
  n <- length(x)
  oldLL <- -2^31 - 1
  for (iter in 1:control$maxiter) {
    cond <- estep(x, Pi, delta, distn, pm, object$pn)
    diff <- cond$LL - oldLL
    if (control$prt) {
      cat("iter =", iter, "\n")
      cat("LL =", formatC(cond$LL, digits=log10(1/1e-5),
                          format="f"), "\n")
      cat("diff =", diff, "\n")
      cat("abs.rel.diff =", abs( diff / oldLL ), "\n\n")
    }
    if (diff < 0 & control$posdiff) stop("Worse log-likelihood on last iteration")
    if (eval(control$converge)) break
    #----  Mstep  ----
    Pi <- diag(1/apply(cond$v, MARGIN = 2, FUN = sum)) %*%
      apply(cond$v, MARGIN = c(2, 3), FUN = sum)
    if (object$nonstat) delta <- cond$u[1, ]
    else delta <- compdelta(Pi)
    pm <- eval(Mstep)
    oldLL <- cond$LL
  }
  object$delta <- delta
  object$Pi <- Pi
  object$u <- cond$u
  object$v <- cond$v
  object$pm <- pm
  object$LL <- cond$LL
  object$iter <- iter
  object$diff <- diff
  return(object)
}
