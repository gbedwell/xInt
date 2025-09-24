#' Mutate Sequences
#'
#' Mutate NGS reads at a defined error rate.
#'
#' @param seq The sequence to be mutated.
#' @param mismatch.prob The mismatch probability. Defaults to 0.0025.
#' @param indel.prob The indel probability. Defaults to 2.5E-5.
#' @param max.bp The maximum read length.
#'
#' @return A mutated DNAString object.
#'
mutate_sequences <- function(seq, mismatch.prob = 0.0025, indel.prob = 2.5E-5, max.bp = 150){
  bases <- strsplit(as.character(seq), "")[[1]]

  mismatch_sites <- which(runif(length(bases)) < mismatch.prob)
  if(length(mismatch_sites) > 0){
    possible_bases <- c("A", "C", "G", "T")
    for (pos in mismatch_sites) {
      bases[pos] <- sample(setdiff(possible_bases, bases[pos]), 1)
    }
  }

  indel_sites <- which(runif(length(bases)) < indel.prob)
  if(length(indel_sites) > 0){
    for (pos in sort(indel_sites, decreasing = TRUE)) {
      if (runif(1) < 0.5) {
        bases <- bases[-pos]
      } else {
        bases <- append(bases, sample(c("A", "C", "G", "T"), 1), after = pos)
      }
    }
  }

  bases <- paste0(bases, collapse = "")
  bases <- substr(bases, 1, max.bp)
  as.character(DNAString(bases))
}
