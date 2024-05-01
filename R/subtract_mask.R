#' Subtract Mask
#'
#' Subtract masked regions of the genome from ranges of interest. Takes GRanges objects as inputs. Will append all reference metadata columns except the revmap column (if present) to the subtracted ranges.
#'
#'@param reference GRanges object corresponding to the regions-of-interest.
#'@param mask GRanges object corresponding to the masked regions.
#'
#'@return A GRanges object containing the subtracted coordinates.
#'
#'@import GenomicRanges
#'@import S4Vectors
#'
#'@export
#'
subtract_mask <- function( reference, mask ){

  diffranges <- GenomicRanges::subtract( x = reference, y = mask,
                                         minoverlap = 1L, ignore.strand = TRUE )
  diffranges <- unlist( diffranges )
  overlaps <- findOverlaps( query = reference, subject = diffranges,
                            minoverlap = 1L, ignore.strand = TRUE )

  n <- names( mcols( reference ) )
  n <- n[ n != "revmap" ]

  for( i in n ){
    mcols( diffranges )[ ,i ] <- NA
    mcols( diffranges )[ ,i ][ subjectHits( overlaps ) ] <-
      mcols( reference )[ ,i ][ queryHits( overlaps ) ]
  }

  return( diffranges)
}
