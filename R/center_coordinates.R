#' Center Coordinates
#'
#' Centers the input site coordinates relative to the defined target site duplication.
#' In the case of odd TSDs, returns a 1 bp interval. For even TSDs, a 2 bp interval is returned.
#'
#' @param sites A SiteList object.
#' @param current.start The position in the target site duplication currently described by start.
#' This is used for centering the site coordinates.
#' @param tsd The total length of the target site duplication.
#' This is used for centering the site coordinates.
#' @param genome.obj The genome object of interest. When not NULL,
#' the validity of the new sites are checked against the genome coordinates.
#' Any centered coordinates that violate genome boundaries are removed.
#'
#' @return A GRanges object containing the centered genomic coordinates.
#'
#' @import methods
#' @import GenomicRanges
#'
center_coordinates <- function(sites, current.start = 1, tsd = 5, genome.obj = NULL){

  if(!validObject(sites)){
    stop("sites is not a valid SiteList object.",
          call. = FALSE)
  }

  ll <- lapply(
    X = sites,
    FUN = function(x){
      center <- (tsd + 1) / 2
      diff <- center - current.start
      a <- floor(diff)
      b <- ceiling(diff)

      plus <- x[strand(x) == "+"]
      end(plus) <- start(plus) + b
      start(plus) <- start(plus) + a

      minus <- x[strand(x) == "-"]
      start(minus) <- end(minus) - b
      end(minus) <- end(minus) - a

      centered <- sort(c(plus, minus), ignore.strand = TRUE)

      if(!is.null(genome.obj)){
        outs <- bound_check(
          fragments = centered,
          genome.obj = genome.obj,
          include.lower = TRUE
          )

        if(length(outs) != 0){
          warning("Centered coordinates " , paste(outs, collapse=", "),
                  " contain out of bounds ranges. These coordinates will be removed.",
                  call. = FALSE)

          centered <- centered[-outs]
          }
        }
        return(centered)
      }
   )

  # ll <- SiteList(ll)

  return(ll)
}
