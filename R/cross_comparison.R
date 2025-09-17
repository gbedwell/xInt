#' Calculate Cross-Sample Overlap
#'
#' Calculates the number of overlapping sites between each sample in a SiteListObject.
#'
#' @param sites A SiteListObject object.
#' @param overlap.coef One of 'overlap' (Szymkiewicz–Simpson overlap coefficient),
#' 'jaccard' (Jaccard index), 'dice', (Dice-Sørensen coefficient), or
#' 'none' (the number of overlapping sites is returned).
#' @param return.sites Returns a GRanges object holding the sites present in 2 or more datasets.
#'
#' @return A matrix of counts denoting the number of overlapping sites between each sample.
#'
#' @examples
#' cross_comparison(sites = sites)
#'
#' @import GenomicRanges
#'
#' @export
#'
cross_comparison = function(sites,
                            overlap.coef = c("overlap", "jaccard", "dice", "none"),
                            return.sites = FALSE){

  if(!validObject(sites)){
    stop("sites is not a valid SiteList object.",
         call. = FALSE )
  }

  uniq.sites <- lapply(sites, unique)

  if(isTRUE(return.sites)){

    all.sites <- unlist(GRangesList(uniq.sites@sites))
    site.counts <- table(all.sites)
    duplicate.sites <- site.counts[which(site.counts >= 2)]
    joined.sites <- paste(seqnames(all.sites), start(all.sites), strand(all.sites), sep=":")
    shared.sites <- unique(all.sites[which(joined.sites %in% names(duplicate.sites))])

    return(shared.sites)

    } else{
    overlap.coef <- match.arg(overlap.coef)

    site.count <- lapply(
      X = seq_along(uniq.sites),
      FUN = function(x){
        tmp.counts <- lapply(
          X = seq_along(uniq.sites),
          FUN = function(y){
            a <- sites[[y]]
            b <- sites[[x]]

            res <- sum(a %in% b)

            if(overlap.coef != "none"){
              if(overlap.coef == "overlap"){
                denom <- min(length(a), length(b))
              } else if(overlap.coef == "jaccard"){
                denom <- length(a) + length(b) - res
              } else if(overlap.coef == "dice"){
                denom <- length(a) + length(b)
                res <- 2 * res
              }
              res <- res / denom
            }
            return(res)
          }
        )
        do.call(c, tmp.counts)
      }
    )
    mat <- do.call(rbind, site.count)
    colnames(mat) <- names(sites)
    rownames(mat) <- names(sites)
    return(mat)
    }
}
