#' Make Sequence Logos
#'
#' Generates sequence logos or consensus matrices from SiteList objects.
#'
#' @param sites A SiteList object.
#' @param seq.len The desired length of the expanded sequences.
#' @param genome.obj The BSgenome object for the genome of interest.
#' @param ignore.strand Boolean. Whether or not to ignore strandedness when extracting sequences.
#' Defaults to FALSE. If TRUE, the sequence from the forward strand is returned.
#' @param current.start The position in the target site duplication currently described by the start coordinates in sites.
#' This is used internally for centering the integration site coordinates.
#' @param tsd The total length of the target site duplication.
#' This is used along with current.start for centering the integration site coordinates.
#' @param n.samp If not NULL, the number of sites to sample from a given dataset. This can significantly speed up processing time,
#' at the expense of including every site in the sequence logo.
#' @param return.plot Boolean.
#' Whether or not to return logos for each dataset in sites instead of the consensus matrix.
#' Defaults to TRUE.
#' @param ... Used to pass additional plotting parameters to ggseqlogo.
#' Examples include ncol and nrow.
#' @param wrap Boolean. Used with return.plot.
#' Whether or not to wrap the returned plots into a single output.
#' When TRUE, arguments ncol and nrow can be used to format the output.
#'
#' @return If return.plot = TRUE, a list of ggplot2 objects containing sequence logos
#' for the mapped integration site datasets in sites.
#' If return.plot = FALSE, a list of consensus matrices for each integration site dataset.
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hs1)
#' data(sites2)
#' make_logo(sites = sites2,
#'           genome.obj = Hsapiens)
#'
#' @import ggplot2
#' @import ggseqlogo
#' @import GenomicRanges
#' @import Biostrings
#'
#' @export
#'
make_logo <- function(sites, seq.len = 24, genome.obj, ignore.strand = FALSE,
                      current.start = 1, tsd = 5, n.samp = NULL, return.plot = TRUE, wrap = FALSE, ... ){

  if(!validObject(sites)){
    stop("sites is not a valid SiteListObject.",
         call. = FALSE)
  }

  if(seq.len %% 2 != 0){
    warning("seq.len cannot be odd. Subtracting 1 from seq.len.")
    seq.len <- seq.len - 1
    }

  expanded.sites <- expand_sites(
    sites = sites@sites,
    seq.len = seq.len,
    genome.obj = genome.obj,
    current.start = current.start,
    tsd = tsd
    )

  if(ignore.strand) {
    expanded.sites <- lapply(
      X = expanded.sites,
      FUN = function(x){
        strand(x) <- "+"
        return(x)
        })}

  if (!is.null(n.samp)) {
    expanded.sites <- lapply(
      X = expanded.sites,
      FUN = function(x){
        if (length(x) > n.samp){
          sample(x = x, size = n.samp, replace = FALSE)
        } else {x}})}

  if(!isTRUE(return.plot)) {
    seqs <- getSeq(
      x = genome.obj,
      names = as(expanded.sites, "GRangesList"),
      as.character = FALSE
      )

    mat.list <- lapply(
      X = seqs,
      FUN = function(x){
        consensusMatrix(x, baseOnly = TRUE)
        })

    return(mat.list)
  } else{
    # For ggseqlogo to apply the small sample correction to
    # the calculated entropy, the sequences themselves must be provided.
    # Providing the PWM itself will not apply any correction.

    seqs <- getSeq(
      x = genome.obj,
      names = as(expanded.sites, "GRangesList"),
      as.character = TRUE 
      )

    center <- (tsd + 1) / 2

    if(tsd %% 2 == 1){
      mp <- (seq.len / 2)
      zp <- mp - (center - 1)
    } else{
      mp <- seq.len / 2
      zp <- mp - (floor(center) - 1)
    }

    p.labs <- seq_len(seq.len) - zp

    if(!isTRUE(wrap)){
      logo.p <- lapply(
        X = seq_along(seqs),
        FUN = function(x){
          ss <- seqs[[x]]

          ggseqlogo( data = ss, seq_type = "dna", ... ) +
            theme_bw() +
            theme(axis.text.y=element_text(size=14),
                  axis.text.x=element_text(size=12),
                  axis.title=element_text(size=16)) +
            scale_x_continuous(breaks=seq( 1, seq.len, 1 ),
                                labels=p.labs,
                                expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            labs(x="Position", title=names(seqs)[x])
          }
        )
    } else{
      logo.p <- ggseqlogo(data = as.list(seqs), seq_type = "dna", ...) +
        theme_bw() +
        theme(axis.text.y=element_text(size=14),
              axis.text.x=element_text(size=12),
              axis.title=element_text(size=16),
              strip.text=element_text(size=12)) +
        scale_x_continuous(breaks=seq( 1, seq.len, 1 ),
                           labels=p.labs,
                           expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        labs(x="Position")
    }

    return(logo.p)
  }
}
