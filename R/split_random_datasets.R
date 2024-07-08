#' Split Random Datasets
#'
#' Splits a list of large random fragment position datasets by chromosome subsets and combines
#' matching subsets from all datasets into a single object.
#'
#'@param p.list A list of datasets holding the positions of random fragments.
#'@param genome.obj The BSgenome object of interest.
#'@param n.cores The number of cores to use.
#'@param n.split The number of chromosome subsets to generate.
#'@param directory.path The path of the output directory.
#'
#'@return Saves to disk compressed RData objects containing random fragment positions for subsets of chromosomes.
#'Note that all saved objects contain a variable of the same name.
#'
#'
#'@import methods
#'
#'@export
#'
split_random_datasets <- function( p.list, genome.obj, n.cores = 1, n.split = 5, directory.path = "." ){

  if( !is.list( p.list ) && !is( p.list, "GRangesList" ) ){
    stop( "p.list must be a list of GRanges objects or a GRangesList.",
          call. = FALSE )
  }

  sl <- seqlengths( genome.obj )
  seqs <- split( names( sl ), ceiling( seq_along( names( sl ) ) / n.split ) )

  ll <- mclapply(
    X = seqs,
    FUN = function(x){
      tmp <- lapply( X = p.list,
                     FUN = function(y){
                       y[seqnames(y) %in% x]
                     }
      )
      tmp <- sort( do.call( c, tmp ), ignore.strand = TRUE )
    },
    mc.cores = n.cores
  )

  invisible(
    lapply(
      X = seq_along(ll),
      FUN = function(x){
        sites <- ll[[x]]
        chrs <- seqlevels(sites)
        chr.from <- chrs[1]
        chr.to <- chrs[length(chrs)]
        save( sites,
              file = paste0( directory.path, "/combined_sites_", chr.from, "to", chr.to, ".RData.gz" ),
              compression_level = 6 )

      }
    )
  )
}
