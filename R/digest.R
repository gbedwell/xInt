#' Simulate Restriction Digestion
#'
#' Simulate restriction digestion of DNA using arbitrary restriction enzymes. DNA to be digested is input as a list of DNAStrings. Restriction sites are defined in a character vector holding the target recognition sequences.
#'
#'@param string.list The output of get_chromosome_seqs() or a comparable list of DNAStrings.
#'@param re.sites A character vector of restriction enzyme recognition sequences.
#'@param cut.after A single integer that defines the nucleotide after which the DNA is cleaved. Defaults to 1.
#'
#'@import Biostrings
#'@import GenomicRanges
#'
#'@returns A list of cut site positions across every chromosome for every defined recognition sequence.
#'
#'@examples
#'\dontrun{
#'digest( string.list = chromosome.sequences, re.sites = c("GGATCC", "TCTAGA") )
#'}
#'
#'@export
#'
digest <- function( string.list, re.sites, cut.after = 1 ){
  seq.obj <- string.list

  for.pattern <- DNAStringSet( re.sites, use.names = TRUE )
  rev.pattern <- complement( DNAStringSet( re.sites, use.names=TRUE ) )

  positions <- lapply(X = seq_along( seq.obj ),
                      FUN = function(x) {
                        for.pos <- matchPDict( for.pattern, seq.obj[[x]] )
                        for.pos <- GRanges( seqnames = names( seq.obj[x] ),
                                                              ranges = unlist( for.pos ),
                                                              strand = "+" )
                        start( for.pos ) <- start( for.pos ) + cut.after
                        end( for.pos ) <- start( for.pos )

                        rev.pos <- matchPDict( rev.pattern, seq.obj[[x]] )
                        rev.pos <- GRanges( seqnames = names( seq.obj[x] ),
                                                              ranges = unlist( rev.pos ),
                                                              strand = "-" )
                        start( rev.pos ) <- end( rev.pos ) - cut.after
                        end( rev.pos ) <- start( rev.pos )

                        pos.list <- do.call( c, list( for.pos, rev.pos ) )

                        return( pos.list ) }
                      )
  return( positions )
}
