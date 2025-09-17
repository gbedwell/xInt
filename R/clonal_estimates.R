#' Clonal Statistics of Integration Site Datasets
#'
#' @param sites A SiteList object
#' @param conditions A character vector of conditions corresponding to the individual datasets.
#' Must be equal to length(sites) and must be in the appropriate order.
#' When not NULL (default), replicates from the same condition are pooled.
#' @param score.col The column in the SiteList object holding site count values.
#' This is included for instances when e.g., site abundance is estimated by MLE methods.
#' When not NULL, this column should be present in the SiteList object.
#'
#' @return A data frame containing the Gini coefficient, Shannon index, Pielou index, and Simpson index for each sample/condition.
#'
#' @importFrom stats p.adjust qnorm pnorm
#'
clonal_estimates <- function(sites, conditions = NULL, score.col = NULL) {
  
  if(!validObject(sites)) {
    stop("sites is not a valid SiteList object.", call. = FALSE)
  }

  if(!is.null(conditions) && length(conditions) != length(sites)) {
    stop("conditions must have the same length as sites.", call. = FALSE)
  }

  if(!is.null(score.col) && !score.col %in% colnames(mcols(sites@sites[[1]]))) {
    stop("score.col not found in sites.", call. = FALSE)
  }

  if(!is.null(conditions)) {
    cond.list <- lapply(
      X = unique(conditions),
      FUN = function(x) {
        which(conditions == x)
      }
    )

    names(cond.list) <- unique(conditions)

    sites <- collapse_sites(sites, cond.list)
    conditions <- unique(conditions)
  } else {
    conditions <- names(sites)
  }
  
  clonality <- lapply(
    X = seq_along(sites),
    FUN = function(i) {
      gr <- sites[[i]]
      nm <- names(sites)[i]
      if(!is.null(score.col)) {
        site.counts <- round(mcols(gr)[[score.col]])
      } else {
        site.counts <- table(gr)
      }

      gini <- calculate_gini(site.counts)
      shannon <- calculate_shannon(site.counts)

      out <- data.frame(
        name = nm,
        gini = gini,
        shannon = shannon[[1]],
        pielou = shannon[[2]],
        simpson = shannon[[3]]
      )
    }
  )

  return(do.call(rbind, clonality))
}

calculate_gini <- function(x) {
  # Sort the values
  x <- sort(as.numeric(x))
  n <- length(x)
  
  # If all values are 0 or there's only one value, return 0
  if (n <= 1 || all(x == 0)) {
    return(0)
  }
  
  cumul.x <- cumsum(x)
  total.x <- cumul.x[n]
  
  # If total is 0, return 0
  if (total.x == 0) {
    return(0)
  }
  
  # Calculate using the formula based on the Lorenz curve
  B <- sum(cumul.x) / (n * total.x)
  gini <- 1 - 2 * B
  
  return(gini)
}

calculate_shannon <- function(x) {
  # If all values are 0 or there's only one value, return 0
  if (length(x) <= 1 || all(x == 0)) {
    return(list(shannon = 0, normalized.shannon = 0))
  }
  
  p <- x / sum(x)
  p <- p[p > 0]  # Remove zeros to avoid log(0); should not happen
  
  # Quantifies both the uniqueness of integration sites (increases with dataset depth) and evenness of their distribution.
  shannon <- -sum(p * log(p))

  # Normalizes shannon (H) between 0 and 1.
  # Measures evenness -- how integration events are distributed across different sites.
  # Lower values indicate dominance of a small number of sites.
  pielou <- shannon / log(length(p))

  # Represents the probability that randomly selected sites are from different positions.
  simpson <- 1 - sum(p^2)      
  
  return(list(shannon = shannon, pielou = pielou, simpson = simpson))
}