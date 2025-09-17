#' Subtract Artifact Sites
#'
#' Subtract potentially artifactual integration sites from all other datasets.
#'
#' @param site.list A list of GRanges objects or a GRangesList containing integration site datasets.
#' @param artifact.datasets A character vector containing the names of artifactual datasets in site.list
#' or a numeric vector containing the indices of artifactual datasets in site.list.
#' These datasets are removed from the output.
#' @param artifact.sites A single GRanges object containing the artifactual sites to be removed from 
#' all entries in site.list.
#'
#' @return A list of GRanges objects or a GRangesList containing the subtracted integration site datasets.
#'
#' @import methods
#' @import GenomicRanges
#'
#' @export
#'
subtract_artifacts <- function(site.list,
                               artifact.datasets = NULL,
                               artifact.sites = NULL){

  if( !validObject( site.list ) ){
    stop("site.list is not a valid SiteListObject.",
         call. = FALSE)
  }

  if(is.null(artifact.datasets) && is.null(artifact.sites)){
    stop("artifact.datasets and artifact.sites cannot both be NULL.", call. = FALSE)
  }

  if(!is.null(artifact.datasets)){
    if(!is.character(artifact.datasets)){
        if(is.numeric(artifact.datasets)){
          artifact.datasets <- names(site.list)[artifact.datasets]
          } else{
            stop("artifact.datasets must be a numeric or character vector.",
                call. = FALSE)
          }
      }
    from.dataset <- c(
      unlist(as(site.list[which(names(site.list) %in% artifact.datasets)]), "GRangesList"),
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

  sub.sites <- lapply(
    X = site.list[-which(names(site.list) %in% artifact.datasets)],
    FUN = function(x){
      ss <- subtract(x = x, y = to.remove, minoverlap = 1L, ignore.strand = FALSE)
      unlist(as(ss, "GRangesList"))
      }
    )

    if(is(site.list, "GRangesList")){
      sub.sites <- as(sub.sites, "GRangesList")
    }

    return(sub.sites)
}
