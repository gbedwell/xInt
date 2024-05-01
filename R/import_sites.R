#' Import Integration Site Coordinates
#'
#' Import a set of integration site datasets from a defined path.
#'
#'@param path The path of the mapped integration site datasets.
#'@param pattern The extension of the integration site dataset files. Defaults to ".bed".
#'@param genome.obj The BSgenome object corresponding to genome of interest.
#'When not NULL, the seqlevels in each IS dataset are compared to the seqlevels of the genome.obj.
#'Seqlevels present in the IS datasets but not present in the genome.obj are removed.
#'@param levels.style The preferred style of seqlevels.
#'@param keep.metadata Boolean. Whether or not the keep metadata contained in the imported datasets. Defaults to FALSE.
#'@param return.GRangesList Boolean. Whether or not the return the list of imported integration site datasets as a GRangesList. Defaults to FALSE.
#'
#'@return A list of GRanges objects or a GRangesList containing mapped integration sites.
#'
#'@import methods
#'@import rtracklayer
#'@import GenomeInfoDb
#'@import GenomicRanges
#'
#'@export
#'
import_sites <- function( path, pattern = ".bed", genome.obj = NULL, levels.style = NULL,
                          keep.metadata = FALSE, return.GRangesList = FALSE ){

  ff <- list.files( path = path,
                    full.names = TRUE,
                    pattern = pattern )

  sites <- lapply( X = ff,
                   FUN = function(x){
                     gr <- import( con = x )

                     if( !is.null( levels.style ) ){
                       seqlevelsStyle( gr ) <- levels.style
                       }

                     if( !keep.metadata ){
                       mcols(gr) <- NULL
                     }

                     return(gr)
                     }
                   )

  nn <- gsub( paste( ".*/(.*?)\\", pattern, sep = "" ), "\\1", ff )

  if( any( duplicated( nn ) ) ){
    stop( "Duplicated filenames found. Make sure that all filenames are unique.",
          "\n",
          "The duplicated filenames are: ", nn[ duplicated( nn ) ], ".",
          call. = FALSE )
  } else{
    names(sites) <- nn
  }

  if( !is.null( genome.obj ) | !missing( genome.obj ) ){
    levs <- lapply( X = sites,
                    FUN = function(x){
                      seqlevels(x)
                    }
    )

    levs <- unique( do.call( c, levs ) )

    out <- setdiff( levs, seqlevels( genome.obj ) )
    levs <- levs[ levs %in% seqlevels( genome.obj ) ]

    if( length( out ) > 0 ){
      warning( "The following seqlevels were removed from the integration site datasets: ",
               paste( out, collapse = ", " ),
               call. = FALSE )
    }

    sites <- lapply( X = sites,
                     FUN = function(x){
                       seqlevels(x) <- levs
                       return(x)
                       }
                     )
    }

  if( return.GRangesList ){
    sites <- as( sites, "GRangesList" )
  }

  return( sites )
}














