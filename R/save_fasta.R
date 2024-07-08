#' Save FASTA
#'
#' Saves simulated sequencing reads as R1 and R2 fasta files.
#'
#'@param reads The output from trim_seqs().
#'@param directory.path The parent directory path.
#'@param prefix The file prefix.
#'@param compress Boolean. Whether or not to compress the output fasta file.
#'@param paired Boolean. Whether or not to save paired-end fasta files.
#'
#'@return None
#'
#'@import Biostrings
#'
save_fasta <- function( reads, directory.path, prefix, paired = TRUE, compress = TRUE ){

  if( isTRUE( compress ) ){
    end <- ".fa.gz"
    } else {
    end <- ".fa"
    }

  R1 <- writeXStringSet( x = reads[[1]],
                         filepath = paste0( directory.path, "/", prefix, "_R1", end ),
                         compress = compress,
                         format = "fasta",
                         width = 1000)

  if( isTRUE( paired ) ){
    R2 <- writeXStringSet( x = reads[[2]],
                           filepath = paste0( directory.path, "/", prefix, "_R2", end ),
                           compress = compress,
                           format = "fasta",
                           width = 1000)
  }
}
