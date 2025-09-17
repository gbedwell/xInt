#' Calculate CPM Values
#'
#' Perform CPM normalization on count values.
#'
#' @param y Input count values
#' @param lib.size A numeric vector of library sizes of each sample.
#' @param log Boolean. Whether or not to return log-CPM. Default TRUE.
#' @param prior.count The constant added to count values to avoid taking log2(0).
#' @param lib.prior The constant added to library size when calculating log-CPM. 
#' When NULL, lib.prior = prior.count.
#'
#' @return A matrix or vector of log-CPM values
#'
calculate_cpm <- function(y, lib.size, log = TRUE, prior.count = 2, lib.prior = NULL) {

  if(!is.numeric(lib.size) && !is.double(lib.size)) {
    stop("lib.size must be double or numeric", call. = FALSE)
  }
  
  if(is.vector(y)) {
    if(length(y) != length(lib.size)) {
      stop("When input is a vector, length(y) must equal length(lib.size).", call. = FALSE)
    }
    y <- matrix(y, nrow = 1)
    return.vector = TRUE
  } else {
    y <- as.matrix(y)
    if(ncol(y) != length(lib.size)) {
      stop("ncol(y) must equal length(lib.size).", call. = FALSE)
    }
    return.vector = FALSE
  }
  
  # Check library sizes
  if (any(is.na(lib.size)) || any(lib.size <= 0)) {
    stop("library sizes should be finite and positive")
  }

  if(is.null(lib.prior)) {
    lib.prior = prior.count
  }
  
  # Calculate CPM
  if(log) {
    # Log CPM with prior count
    cpm <- log2((y + prior.count) * 1e6 / (lib.size + lib.prior))
  } else {
    # Raw CPM
    cpm <- t(t(y) * 1e6 / lib.size)
  }
  
  # Preserve dimension names
  dimnames(cpm) <- dimnames(y)

  if(return.vector) {
    return(as.vector(cpm))
  } else {
    return(cpm)
  }
}