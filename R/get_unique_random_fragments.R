#' Get Unique Random Fragments
#'
#' Extracts uniquely mapped random fragments from bam files generated using simulated random data. Iterates through the bam file in yield.size chunks to save memory. Outputs the random fragment positions in a GRanges object and saves the object to disk as an RData file.
#'
#'@param bam.file The bam file from the alignment
#'@param mapq.cutoff The desired mapq cutoff value. If no quality fltering is desired, set to 0.
#'@param yield.size The number of reads from the bam file to process at a time.
#'@param sort Logical. Whether or not to sort the generated GRanges object. Defaults to FALSE. Ignores strandedness.
#'@param prefix The desired output file prefix. Used to name the output RData file.
#'@param write Boolean. Whether or not to save the output to disk.
#'
#'@return A GRanges object containing the uniquely mapped random fragments.
#'
#'@import Rsamtools
#'@import GenomicAlignments
#'@import GenomicFiles
#'@importFrom dplyr filter
#'@importFrom dplyr mutate
#'@importFrom dplyr select
#'
#'@export
#'
get_unique_random_fragments <- function( bam.file,
                                         mapq.cutoff,
                                         yield.size = NA,
                                         sort = FALSE,
                                         prefix,
                                         write = FALSE ){

  bam <- BamFile( bam.file, asMates = TRUE, yieldSize = yield.size )

  params <- ScanBamParam( flag = scanBamFlag( isProperPair = TRUE,
                                              isDuplicate = FALSE,
                                              isSecondaryAlignment = FALSE ),
                          mapqFilter = mapq.cutoff )

  yield <- function(x){
    readGAlignmentPairs( x, param = params )
  }

  map <- identity
  reduce <- c

  out <- reduceByYield( X = bam, YIELD = yield, MAP = map, REDUCE = reduce )

  df <- data.frame( seqnames( out ), ranges( out ), strand( out ) )
  colnames( df ) <- c( "seqnames", "start", "end", "width", "strand" )

  gr <- GRanges( df )

  if ( isTRUE( sort ) ){
    gr <- sort( gr, ignore.strand = TRUE )
  }

  if ( !isTRUE( write ) ){
    return( gr )
    } else{
    if ( missing( prefix ) ){
      stop( "prefix must be defined out write out the .RData file.",
            call. = FALSE )
      }
      save( gr,
            file = paste0( prefix, "_unique_random_fragments.RData.gz") ,
            compression_level = 6 )
    }
  }














