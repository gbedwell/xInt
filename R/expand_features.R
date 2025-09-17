#' Expand Feature Coordinates
#'
#' Expands feature coordinates by a defined amount.
#'
#'@param features The GRanges object holding the features-of-interest
#'@param upstream The amount to expand (in bp) upstream of the FOI.
#' Features with strand = '*' are treated like strand = '+'.
#'@param downstream The amount to expand (in bp) downstream of the FOI.
#'@param genome.obj The genome object of interest. This is used for checking the validity of
#'the expanded coordinates. Can be NULL.
#'@param ignore.strand Boolean. Whether or not to ignore strandedness. In this case, all features are treated like strand = '+'.
#' Defaults to FALSE
#'
#' @return A GRanges object containing the expanded coordinates.
#'
#' @import methods
#' @import GenomicRanges
#'
#' @export
#' 

expand_features <- function(features, upstream = 0, downstream = 0, genome.obj = NULL, ignore.strand = FALSE) {
  upstream <- abs(upstream)
  downstream <- abs(downstream)  # Add this line

  if(isTRUE(ignore.strand)) {
    og.strand <- strand(features)
    strand(features) <- "*"
  }

  # Fix the syntax error here
  strand.plus <- as.logical(strand(features) == "+" | strand(features) == "*")
  strand.minus <- as.logical(strand(features) == "-")
  
  if(any(strand.plus)) {
    # For + strand: upstream expands start, downstream expands end
    start(features[strand.plus]) <- pmax(1, start(features[strand.plus]) - upstream)
    end(features[strand.plus]) <- end(features[strand.plus]) + downstream
  }
  
  if(any(strand.minus)) {
    # For - strand: upstream expands end, downstream expands start
    start(features[strand.minus]) <- pmax(1, start(features[strand.minus]) - downstream)
    end(features[strand.minus]) <- end(features[strand.minus]) + upstream
  }

  # Fix the genome length checking
  if(!is.null(genome.obj)) {
    seq.lens <- seqlengths(genome.obj)[as.character(seqnames(features))]
    extreme <- which(end(features) > seq.lens)  # Fix: was missing parentheses
    if(length(extreme) > 0) {
      end(features)[extreme] <- seq.lens[extreme]
    }
  } else if(!any(is.na(genome(features)))) {
    features <- trim(features)
  }

  if(isTRUE(ignore.strand)) {
    strand(features) <- og.strand
  }

  return(features)
}