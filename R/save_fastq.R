#' Save FASTQ
#'
#' Saves simulated sequencing reads as R1 and R2 FASTQ files.
#'
#' @param reads The output from trim_seqs().
#' @param directory.path The parent directory path.
#' @param prefix The file prefix.
#' @param compress Boolean. Whether or not to compress the output fasta file.
#' @param paired Boolean. Whether or not to save paired-end fasta files.
#'
#' @return None
#'
#' @import Biostrings
#'
save_fastq <- function(reads, directory.path, prefix, paired = TRUE, compress = TRUE) {

  r1_file <- file.path(directory.path,
                       paste0(prefix, "_R1", ifelse(compress, ".fq.gz", ".fq")))
  r2_file <- file.path(directory.path,
                       paste0(prefix, "_R2", ifelse(compress, ".fq.gz", ".fq")))

  write_fastq <- function(seqs, quals, file) {
    con <- if(compress) gzfile(file, "w") else file(file, "w")
    for(i in seq_along(seqs)) {
      writeLines(c(
        paste0("@", names(seqs)[i]),
        seqs[i],
        "+",
        quals[i]
      ), con)
    }
    close(con)
  }

  write_fastq(reads$R1$seq, reads$R1$qual, r1_file)
  if(isTRUE(paired)){
    write_fastq(reads$R2$seq, reads$R2$qual, r2_file)
  }
}
