#' Simulate Restriction Digestion
#'
#' Simulate restriction digestion of DNA using arbitrary restriction enzymes.
#' DNA to be digested is input as a list of DNAStrings.
#' Restriction sites are defined in a character vector holding the target recognition sequences.
#'
#' @param chr.seqs The output of get_chromosome_seqs() or a comparable list of DNAStrings.
#' @param re.sites A character vector of restriction enzyme recognition sequences.
#' @param cut.after A single integer that defines the nucleotide after which the DNA is cleaved. Defaults to 1.
#'
#'
#' @import Biostrings
#' @import GenomicRanges
#'
#' @return A list of cut site positions across every chromosome for every defined recognition sequence.
#'
digest <- function(chr.seqs, re.sites, cut.after = 1){

  for.pattern <- DNAStringSet(re.sites, use.names = TRUE)
  rev.pattern <- complement(DNAStringSet(re.sites, use.names=TRUE))

  positions <- lapply(X = seq_along(chr.seqs),
                      FUN = function(x) {
                        for.pos <- matchPDict(for.pattern, chr.seqs[[x]])

                        if(sum(do.call(c, lapply(X = for.pos, FUN = function(x){ length(x) }))) == 0){
                          for.pos <- GRanges(seqnames = NULL, ranges = NULL, strand = NULL)
                        } else{
                          for.pos <- GRanges(seqnames = names(chr.seqs[x]),
                                             ranges = unlist(for.pos),
                                             strand = "+")

                          start(for.pos) <- start(for.pos)
                          end(for.pos) <- start(for.pos) + (cut.after - 1)
                        }

                        rev.pos <- matchPDict(rev.pattern, complement(chr.seqs[[x]]))

                        if(sum(do.call(c, lapply(X = rev.pos, FUN = function(x){length(x)}))) == 0){
                          rev.pos <- GRanges(seqnames = NULL, ranges = NULL, strand = NULL)
                        } else{
                          rev.pos <- GRanges(seqnames = names(chr.seqs[x]),
                                             ranges = unlist(rev.pos),
                                             strand = "-")

                          start(rev.pos) <- end(rev.pos) - (cut.after - 1)
                          end(rev.pos) <- end(rev.pos)
                        }

                        pos.list <- unique(unlist(as(list(for.pos, rev.pos), "GRangesList")))
                        pos.list <- sort(pos.list, ignore.strand = TRUE)
                        return(pos.list)
                        }
                     )
  return(positions)
}
