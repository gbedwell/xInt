#' Trim Fragments
#'
#' Trim the generated random fragments to simulate sequencing reads.
#' Creates paired end fragments of a defined maximum read size and with a defined maximum inner distance.
#'
#' @param fragments The DNAStringSet object containing the random fragment sequences.
#' @param min.width The minimum acceptable fragment length (insert size). Defaults to 25.
#' @param max.width The maximum acceptable fragment length (insert size). Defaults to 1200.
#' @param max.bp The maximum read length. Defaults to 150 bp.
#' @param ltr.seq The LTR sequence sequenced in the experiment, given in the 5' to 3' orientation.
#' Should be a character string. Use 'X' for sample-specific barcodes, 'N' for unique per-sequence barcodes.
#' @param linker.seq The linker sequence sequenced in the experiment, given in the 5' to 3' orientation.
#' Should be a character string. Use 'X' for sample-specific barcodes, 'N' for unique per-sequence barcodes.
#' @param mismatch.prob The mismatch rate. Defaults to 0.0025.
#' @param indel.prob The indel rate. Defaults to 2.5E-5.
#' @param offset The ASCII offset of the quality scores. Defaults to 33.
#' @param base.score The quality score value to assign to all bases. Offset is added to this value. Default 40.
#' @param noisy.scores Boolean. Whether or not to introduce simulated noise into quality scores. Defaults to FALSE.
#' @param read.name.prefix The naming prefix for read names. Defaults to "fragment_".
#' @param name.idx.start The starting index for the naming scheme. Defaults to 0.
#' @param use.frag.names Boolean. Whether or not use names from the input fragments for sequence names. Defaults to TRUE.
#' @param paired Boolean. Whether to generate paired-end reads (TRUE) or single-end reads (FALSE). Defaults to TRUE.
#'
#' @return A list containing simulated R1 and R2 reads.
#'
#' @import GenomicRanges
#' @import Biostrings
#'
trim_seqs <- function(fragments, min.width = 25, max.width = 1200, max.bp = 150,
                      ltr.seq, linker.seq, mismatch.prob = 0.0025, indel.prob = 2.5E-5,
                      base.score = 40, offset = 33, as.fasta = FALSE, noisy.scores = FALSE,
                      use.frag.names = TRUE, read.name.prefix = "fragment_", name.idx.start = 0,
                      paired = TRUE){

  if(missing(ltr.seq) || missing(linker.seq)){
    stop("Both ltr.seq and linker.seq must be defined. Define as NULL if no
         LTR or linker sequences should be simulated.")
  }

  if((is.null(ltr.seq) && !is.null(linker.seq)) || (!is.null(ltr.seq) && is.null(linker.seq))){
    stop("ltr.seq and linker.seq must both be NULL or not NULL.")
  }

  filtered <- fragments[width(fragments) >= min.width & width(fragments) <= max.width]

  if(length(filtered) == 0) {
    warning("No fragments passed the filtering criteria.")
    return(list(R1 = NULL, R2 = NULL))
  }

  if(!is.null(ltr.seq) && !is.null(linker.seq)){

    ltr.has.bc <- grepl("X", ltr.seq)
    linker.has.bc <- grepl("X", linker.seq)

    if(ltr.has.bc || linker.has.bc){
      ltr.seq <- make_sample_barcodes(ltr.seq)
      linker.seq <- make_sample_barcodes(linker.seq)
    }

    ltr.has.umi <- grepl("N", ltr.seq)
    linker.has.umi <- grepl("N", linker.seq)

    if(ltr.has.umi || linker.has.umi){
      all.seqs <- make_umis(
        ltr.seq = ltr.seq,
        linker.seq = linker.seq,
        seqs = filtered,
        ltr.has.umi = ltr.has.umi,
        linker.has.umi = linker.has.umi
      )
      appended <- DNAStringSet(all.seqs)
    } else {
      linker.rc <- as.character(reverseComplement(DNAString(linker.seq)))
      appended <- DNAStringSet(paste0(ltr.seq, as.character(filtered), linker.rc))
    }
  } else{
    appended <- filtered
  }

  left <- subseq(
    x = appended,
    start = 1,
    width = pmin(width(appended), max.bp)
    )

  if(paired) {
    right <- subseq(
      x = appended,
      start = pmax(1, width(appended) - max.bp + 1),
      width = pmin(width(appended), max.bp)
      )

    right <- reverseComplement(right)
  }

  if(mismatch.prob > 0 && indel.prob > 0){
    r1.reads <- DNAStringSet(
      vapply(
        X = left,
        FUN = mutate_sequences,
        FUN.VALUE = character(1),
        mismatch.prob = mismatch.prob,
        indel.prob = indel.prob,
        max.bp = max.bp
        )
      )
    if(paired) {
      r2.reads <- DNAStringSet(
        vapply(
          X = right,
          FUN = mutate_sequences,
          FUN.VALUE = character(1),
          mismatch.prob = mismatch.prob,
          indel.prob = indel.prob,
          max.bp = max.bp
          )
        )
      }
    } else{
    r1.reads <- left
    if(paired) {
      r2.reads <- right
    }
    }

  if(use.frag.names){
    names(r1.reads) <- names(filtered)
    if(paired) {
      names(r2.reads) <- names(filtered)
    }
  } else{
    names(r1.reads) <- paste0(read.name.prefix, name.idx.start + seq_along(r1.reads))
    if(paired) {
      names(r2.reads) <- paste0(read.name.prefix, name.idx.start + seq_along(r2.reads))
    }
  }

  if(!isTRUE(as.fasta)){
    r1.widths <- width(r1.reads)
    if(paired) {
      r2.widths <- width(r2.reads)
    }

    if(!noisy.scores){
      scores <- rep(rawToChar(as.raw(base.score + offset)), max.bp)

      r1.scores <- vapply(
        X = r1.widths,
        FUN = function(x){
          paste(scores[1:x], collapse = "")
        },
        FUN.VALUE = character(1)
      )

      if(paired) {
        r2.scores <- vapply(
          X = r2.widths,
          FUN = function(x){
            paste(scores[1:x], collapse = "")
          },
          FUN.VALUE = character(1)
        )
      }
    } else{
      n.seqs <- length(r1.widths)

      r1.scores <- generate_quality_scores(
        n.seqs = n.seqs,
        max.bp = max.bp,
        base.score = base.score
      )

      r1.scores <- substr(r1.scores, 1, r1.widths)

      if(paired) {
        r2.scores <- generate_quality_scores(
          n.seqs = n.seqs,
          max.bp = max.bp,
          base.score = base.score
        )

        r2.scores <- substr(r2.scores, 1, r2.widths)
      }
    }

    if(paired) {
      return(
        list(
          R1 = list(seq = as.character(r1.reads), qual = r1.scores),
          R2 = list(seq = as.character(r2.reads), qual = r2.scores)
        ))
    } else {
      return(
        list(
          R1 = list(seq = as.character(r1.reads), qual = r1.scores)
        ))
    }
    } else {
    if(paired) {
      return(
        list(
          R1 = r1.reads,
          R2 = r2.reads
        ))
    } else {
      return(
        list(
          R1 = r1.reads
        ))
    }
    }
}

make_sample_barcodes <- function(seq){
  seq.chars <- strsplit(seq, "")[[1]]
  x.positions <- which(seq.chars == "X")
  sample.barcode <- sample(c("A", "T", "G", "C"), length(x.positions), replace = TRUE)

  seq.chars[x.positions] <- sample.barcode

  return(seq = paste0(seq.chars, collapse = ""))
}

make_umis <- function(ltr.seq, linker.seq, seqs, ltr.has.umi, linker.has.umi){
  n.seqs <- length(seqs)

  ltr.chars <- strsplit(ltr.seq, "")[[1]]
  ltr.n.positions <- which(ltr.chars == "N")
  linker.chars <- strsplit(linker.seq, "")[[1]]
  linker.n.positions <- which(linker.chars == "N")

  if(ltr.has.umi) {
    ltr.umi.matrix <- sample(c("A", "T", "G", "C"), length(ltr.n.positions) * n.seqs, replace = TRUE)
    dim(ltr.umi.matrix) <- c(n.seqs, length(ltr.n.positions))

    ltr.template <- matrix(rep(ltr.chars, n.seqs), nrow = n.seqs, byrow = TRUE)

    for(i in seq_along(ltr.n.positions)) {
      ltr.template[, ltr.n.positions[i]] <- ltr.umi.matrix[, i]
    }

    ltr.strings <- apply(ltr.template, 1, paste0, collapse = "")
  } else {
    ltr.strings <- rep(ltr.seq, n.seqs)
  }

  if(linker.has.umi) {
    rc_map <- c(A="T", T="A", G="C", C="G", N="N")

    linker.umi.matrix <- sample(c("A", "T", "G", "C"), length(linker.n.positions) * n.seqs, replace = TRUE)
    dim(linker.umi.matrix) <- c(n.seqs, length(linker.n.positions))

    linker.template <- matrix(rep(linker.chars, n.seqs), nrow = n.seqs, byrow = TRUE)

    for(i in seq_along(linker.n.positions)) {
      linker.template[, linker.n.positions[i]] <- linker.umi.matrix[, i]
    }

    linker.strings <- apply(linker.template, 1, paste0, collapse = "")

    linker.rc.strings <- character(n.seqs)
    for(i in 1:n.seqs) {
      rev_chars <- rev(strsplit(linker.strings[i], "")[[1]])
      comp_chars <- rc_map[rev_chars]
      linker.rc.strings[i] <- paste0(comp_chars, collapse = "")
    }
  } else {
    current.linker.str <- paste0(linker.chars, collapse = "")
    linker.rc <- as.character(reverseComplement(DNAString(current.linker.str)))
    linker.rc.strings <- rep(linker.rc, n.seqs)
  }

  all.seqs <- paste0(ltr.strings, as.character(seqs), linker.rc.strings)

  return(all.seqs)
}

generate_quality_scores <- function(n.seqs = 1000,
                                    max.bp = 150,
                                    base.score = 40,
                                    taper.intensity = 0.2,
                                    as.probs = FALSE,
                                    as.ASCII = TRUE){

  if(as.probs && as.ASCII){
    stop("Both as.probs and as.ASCII cannot be TRUE.",
         call. = FALSE)
  }

  template.scores <- rep(as.integer(base.score), max.bp)

  position.effect <- (1 - exp(-(seq(0, 3, length.out = max.bp))))
  position.effect <- position.effect / max(position.effect)

  quality.decrease <- taper.intensity * base.score * position.effect
  template.scores <- as.integer(round(base.score - quality.decrease))

  noise.scale <- seq(1, 2, length.out = max.bp)

  if(as.ASCII){
    ascii_table <- sapply(as.raw(0:40 + 33), rawToChar)
    names(ascii_table) <- as.character(0:40)
  }

  template_matrix <- matrix(rep(template.scores, each = n.seqs), ncol = max.bp)

  noise_matrix <- matrix(rnorm(n.seqs * max.bp, 0, 0.4), ncol = max.bp)
  for(i in 1:max.bp) {
    noise_matrix[,i] <- noise_matrix[,i] * noise.scale[i]
  }

  noisy.scores <- template_matrix + noise_matrix

  if(!as.probs && !as.ASCII) {
    noisy.scores <- matrix(as.integer(round(noisy.scores)), nrow = n.seqs, ncol = max.bp)
  }

  if(as.ASCII){
    quality.strings <- character(n.seqs)
    for(i in 1:n.seqs) {
      rounded <- pmax(20, pmin(40, round(noisy.scores[i,])))
      chars <- ascii_table[as.character(rounded)]
      quality.strings[i] <- paste(chars, collapse = "")
    }
    return(quality.strings)
  } else if(as.probs){
    probs <- 10^(-noisy.scores / 10)
    return(probs)
  } else{
    return(noisy.scores)
  }
}
