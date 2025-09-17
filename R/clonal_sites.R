#' Define Clonal Sites
#'
#' @param sites A SiteList object
#' @param conditions A character vector of conditions corresponding to the individual datasets.
#' Must be equal to length(sites) and must be in the appropriate order.
#' When not NULL (default), replicates from the same condition are pooled.
#' @param min.threshold The minimum threshold for a site to be considered clonal, expressed as a fraction of all sites 
#' or as a fraction of the total estimated abundance. Defaults to 0 (all sites are returned).
#' @param score.col The column in the SiteList object holding site count values.
#' This is included for instances when e.g., site abundance is estimated by MLE methods.
#' When not NULL, this column should be present in the SiteList object.
#'
#' @return A GRanges object
#'
#' @import GenomicRanges
#'
clonal_sites <- function(sites, conditions, min.threshold = 0, score.col = NULL) {
    
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

  clonal.fracs <- lapply(
    X = seq_along(sites),
    FUN = function(i) {
      gr <- sites[[i]]
      nm <- names(sites)[i]
      if(!is.null(score.col)) {
        site.counts <- round(mcols(gr)[[score.col]])
        site.fracs <- site.counts / sum(site.counts)
        gr$fraction <- site.fracs
      } else {
        site.counts <- tabulate(match(gr, unique(gr)))
        site.fracs <- site.counts / sum(site.counts)
        gr <- unique(gr)
        gr$fraction <- site.fracs
      }
      
      gr <- gr[order(gr$fraction, decreasing = TRUE)]
      gr <- gr[which(gr$fraction >= min.threshold)]
      return(gr)
    }
  )

  return(clonal.fracs)
}