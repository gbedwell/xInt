#' Get Random Coordinates
#'
#' Extracts uniquely mapped random fragments from bam files generated using simulated random data.
#' Iterates through the bam file in yield.size chunks to save memory.
#' Outputs the random fragment positions in a GRanges object and saves the object to disk as an RData file.
#'
#'@param bam.file The bam file from the alignment
#'@param paired Boolean. Whether or not the aligned data came from single- or paired-end reads. Defaults to TRUE.
#'@param return.reads Boolean. Whether or not to return the start/end coordinates for each read or for each mapped fragment.
#'@param mapq.cutoff The desired mapq cutoff value.
#'@param simple.cigar Boolean. Identical to simpleCigar in Rsamtools::ScanBamParam().
#'@param tag A vector of tags on which to filter the input BAM file. To omit, set tag=character(0).
#'@param tagFilter List of possible values for the given tags. To omit entirely, set tagFilter=list().
#'@param yield.size The number of reads from the bam file to process at a time.
#'@param sort Boolean. Whether or not to sort the generated GRanges object. Defaults to FALSE. Ignores strandedness.
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
get_random_coordinates <- function( bam.file,
                                paired,
                                return.reads = FALSE,
                                mapq.cutoff = 0,
                                simple.cigar = FALSE,
                                tag = "NM",
                                tagFilter = list(NM = 0),
                                yield.size = NA,
                                sort = FALSE,
                                prefix,
                                write = FALSE ){

  bam <- BamFile( bam.file, asMates = paired, yieldSize = yield.size )

  if( !isTRUE( return.reads ) ){
    if( isTRUE( paired ) ){
      params <- ScanBamParam( flag = scanBamFlag( isProperPair = TRUE,
                                                  isDuplicate = FALSE,
                                                  isSecondaryAlignment = FALSE,
                                                  isUnmappedQuery = FALSE ),
                              mapqFilter = mapq.cutoff,
                              simpleCigar = simple.cigar,
                              tag = tag,
                              tagFilter = tagFilter )

      yield <- function(x){ readGAlignmentPairs( x, param = params ) }
    } else{
      warning( "Cannot return fragment coordinates for unpaired data.",
               "\n",
               "Treating return.reads as TRUE.",
               call. = FALSE )
    }
  } else{
    params <- ScanBamParam( flag = scanBamFlag( isDuplicate = FALSE,
                                                isSecondaryAlignment = FALSE,
                                                isUnmappedQuery = FALSE ),
                            mapqFilter = mapq.cutoff,
                            simpleCigar = simple.cigar,
                            tag = tag,
                            tagFilter = tagFilter )
    yield <- function(x){
      readGAlignments( x, param = params )
    }
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
            file = paste0( prefix, "_unique_fragments.RData.gz") ,
            compression_level = 6 )
      }
  }














