#' Pool Replicate Datasets
#'
#' Given a SiteList object containing replicates from various conditions,
#' pool the replicates into a single, larger dataset while keeping different conditions separate.
#'
#' @param sites A SiteList object.
#' @param group.index.list A list of index vectors enumerating which sites elements correspond to which conditions.
#' Each group.index.list element should correspond to a unique condition.
#' @param sorted Boolean. Whether or not to sort the concatenated data. 
#' Ignores strand. Defaults to TRUE.
#'
#' @return A list of GRanges objects or a GRangesList of the combined coordinates.
#'
#' @examples
#' data(sites)
#' collapse_sites(sites = sites,
#'                group.index.list = list(A = 1:4,
#'                                        B = 5:9,
#'                                        C = 10)
#'
#' @import methods
#' @import GenomicRanges
#'
#' @export
#'
collapse_sites <- function(sites, group.index.list, sorted = TRUE){

  if(!is.list(group.index.list)){
    stop("Group index values must be given as a list.",
          call. = FALSE)
    }

  if(!validObject(sites)){
    stop("sites is not a valid SiteList object.",
         call. = FALSE)
    }

  ll <- lapply(
    X = group.index.list,
    FUN = function(x){
      ind <- x
      out <- unlist(GRangesList(sites@sites[ind]))
      if(sorted) {
        out <- sort(out, ignore.strand = TRUE)
      }
      return(out)
    }
  )

  if(is.null(names(ll))){
    warning("No group names were provided in
            'group.index.list. The returned list
            will therefore not be named.",
            call. = FALSE)
    }

  ll <- SiteList(ll)

  return(ll)
}
