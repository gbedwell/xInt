#' Generate Random Sites
#'
#' Generates randomly positioned integration sites within the genome of interest.
#'
#' @param n.sites The number of sites to generate.
#' @param genome.obj The BSgenome object of interest.
#' @param n.clones The number of clonal sites to generate. Defaults to 0.
#' @param clone.fracs A numeric vector defining the fraction of sites occupied by each clone.
#' Defaults to 0.
#' @param chr.omit A character vector of chromosomes to omit. Defaults to c("chrM", "MT")
#' Must be of length == n.clones, all values must be between 0 and 1, and the sum of values must be <= 1.
#' @param clone.chroms A vector denoting the chromosome harboring each clonal site. Defaults to NULL.
#' If not NULL, must be of length equal to n.clones.
#' @param clone.positions A vector denoting the position of each clonal site. Defaults to NULL.
#' If not NULL, must be of length equal to n.clones.
#' @param clone.strands A vector denoting the strand of each clonal site. Defaults to NULL.
#' If not NULL, must be of length equal to n.clones.
#'
#' @return A GRanges object containing the integration site positions.
#'
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#'
random_sites <- function(n.sites, genome.obj, n.clones = 0, clone.fracs = 0, chr.omit = c("chrM", "MT"),
                         clone.chroms = NULL, clone.positions = NULL, clone.strands = NULL){

  if(n.clones != 0){
    if(length(clone.fracs) != n.clones){
      stop("length(clone.fracs) != n.clones.",
           call. = FALSE)
    }
    if(any(clone.fracs <= 0) || any(clone.fracs > 1)){
      stop("All clone.fracs values must be > 0 and <= 1.",
           call. = FALSE)
    }
    if(sum(clone.fracs) > 1 || sum(clone.fracs) == 0){
      stop("sum(clone.fracs) must be <= 1 and > 0.",
           call. = FALSE)
    }
    # TO-DO: add checks for clone.starts, clone.chroms, and clone.strands
  }

  chr.names <- seqlevels(genome.obj)[!seqlevels(genome.obj) %in% chr.omit]
  chr.lengths <- seqlengths(genome.obj)[!seqlevels(genome.obj) %in% chr.omit]
  chr.weights <- chr.lengths / sum(chr.lengths)

  n.clone <- floor(clone.fracs * n.sites)
  n.remaining <- n.sites - sum(n.clone)

  chr.vec <- sample(x = chr.names, size = n.remaining, prob = chr.weights, replace = TRUE)
  site.pos <- vapply(X = chr.vec,
                     FUN = function(x){
                       cl <- chr.lengths[[x]]
                       sample(x = seq_len(cl), size = 1)
                       },
                     FUN.VALUE = numeric(1)
  )
  site.strand <- sample(x = c("+", "-"), size = n.remaining, replace = TRUE)
  sites.gr <- GRanges(
    seqnames = chr.vec,
    ranges = IRanges(start = site.pos, width = 1),
    strand = site.strand
    )

  if(n.clones > 0){
    clone.sites <- lapply(
      X = seq_along(n.clone),
      FUN = function(x){
        num.sites <- n.clone[x]
        if(is.null(clone.chroms)){
          clone.chr <- rep(sample(x = chr.names, size = 1, prob = chr.weights), num.sites)
        } else{
          clone.chr <- rep(clone.chroms[x], num.sites)
        }
        cl <- chr.lengths[clone.chr[1]]
        if(is.null(clone.positions)){
          clone.pos <- rep(sample(x = seq_len(cl), size = 1), num.sites)
        } else{
          clone.pos <- rep(clone.positions[x], num.sites)
        }
        if(is.null(clone.strands)){
          clone.strand <- rep(sample(x = c("+", "-"), size = 1), num.sites)
        } else{
          clone.strand <- rep(clone.strands[x], num.sites)
        }

        gr <- GRanges(
          seqnames = clone.chr,
          ranges = IRanges(start = clone.pos, width = 1),
          strand = clone.strand
        )
      }
    )

    # Suppress warnings that no seqlevels are in common
    suppressWarnings(
      clone.sites <- do.call(c, clone.sites)
    )

    suppressWarnings(
      sites.gr <- c(sites.gr, clone.sites)
    )
  }

  return(sites.gr)

}
