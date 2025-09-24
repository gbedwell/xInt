#' Subtract Artifact Sites
#'
#' Subtract potentially artifactual integration sites from all other datasets.
#'
#' @param sites A SiteList object.
#' @param artifact.datasets A character vector containing the names of artifactual datasets in sites
#' or a numeric vector containing the indices of artifactual datasets in sites.
#' These datasets are removed from the output.
#' @param artifact.sites A single GRanges object containing the artifactual sites to be removed from 
#' all entries in sites.
#'
#' @return A list of GRanges objects or a GRangesList containing the subtracted integration site datasets.
#'
#' @import methods
#' @import GenomicRanges
#'
#' @export
#'
subtract_artifacts <- function(sites,
                               artifact.datasets = NULL,
                               artifact.sites = NULL){

  if(!validObject(sites) ){
    stop("sites is not a valid SiteList object.",
         call. = FALSE)
  }

  if(is.null(artifact.datasets) && is.null(artifact.sites)){
    stop("artifact.datasets and artifact.sites cannot both be NULL.", call. = FALSE)
  }

  if(!is.null(artifact.datasets)){
    if(!is.character(artifact.datasets)){
        if(is.numeric(artifact.datasets)){
          artifact.datasets <- names(sites)[artifact.datasets]
          } else{
            stop("artifact.datasets must be a numeric or character vector.",
                 call. = FALSE)
          }
      }
    from.dataset <- c(
      unlist(as(sites[which(names(sites) %in% artifact.datasets)]), "GRangesList"),
    )
  } else{
    from.dataset <- GRanges()
  }
  

  if(is.null(artifact.sites)){
    from.gr <- GRanges()
  } else{
    from.gr <- artifact.sites
  }
  
  to.remove <- c(from.dataset, from.gr)

  if(length(to.remove) == 0) {
    warning("No artifact sites to remove.", call. = FALSE)
    return(sites)
  }

  if(!is.null(artifact.datasets)) {
    sites <- sites[-which(names(sites) %in% artifact.datasets)]
  }

  sub.sites <- lapply(
    X = sites,
    FUN = function(x){
      ss <- subtract(x = x, y = to.remove, minoverlap = 1L, ignore.strand = FALSE)
      unlist(as(ss, "GRangesList"))
      }
    )
  
  sub.sites <- sub.sites[lengths(sub.sites) > 0]

  if(is(sites, "GRangesList")){
    sub.sites <- as(sub.sites, "GRangesList")
  }

  return(sub.sites)
}
