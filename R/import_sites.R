#' Import Integration Site Coordinates
#'
#' Import a set of integration site datasets from a defined path.
#'
#' @param path The path to the mapped integration site datasets.
#' Ideally, the path should only contain the IS datasets of interest.
#' @param files A character vector of file names to import.
#' Allows flexibility if not all datasets are in the same path or
#' if not all files in a given directory are of interest.
#' Should not be used if 'path' is defined.
#' @param pattern The extension of the integration site dataset files. Defaults to ".bed".
#' @param ric Path to a directory or file containing random integration control sites.
#' @param genome.obj The BSgenome object corresponding to genome of interest.
#' When not NULL, the seqlevels in each IS dataset are compared to the seqlevels of the genome.obj.
#' Seqlevels present in the IS datasets but not present in the genome.obj are removed.
#' @param levels.style The preferred style of seqlevels.
#' @param keep.metadata Boolean.
#' Whether or not the keep metadata contained in the imported datasets.
#' Defaults to FALSE.
#' @param sort Boolean. Whether or not to sort the imported data.
#' Sorts seqlevels first, then ranges. Ignores strand.
#' @param rm.string A common string/pattern to remove from all imported filenames.
#' Files containing this string will still get imported, but the string itself will be removed.
#' Defaults to "_sites". Set to NULL to ignore this argument entirely.
#' @param ignore.chr Chromosome names to remove from the imported sites.
#' Defaults to "c("chrM", "MT")".
#' @param ignore.string Filenames containing ignore.string will not be imported. 
#' This allows greater selectivity in which files get imported.
#' Multiple strings can be defined: "stringA|stringB".
#' Defaults to "unique".
#' @param expand.by.count A character vector of length 1 defining the column holding the estimated site abundance values.
#' When not NULL (default), the sites in each file are expanded by the corresponding value in the defined column.
#' @param as.fragments Boolean. If TRUE, returns a list of GRanges objects instead of a SiteList object.
#' Defaults to FALSE.
#' 
#'
#' @return A SiteList object containing mapped integration sites or a list of GRanges objects.
#'
#' @import methods
#' @import rtracklayer
#' @import GenomeInfoDb
#' @import GenomicRanges
#'
#' @export
#'
import_sites <- function(path, files, pattern = ".bed", genome.obj = NULL, levels.style = NULL,
                         keep.metadata = FALSE, sort = TRUE, rm.string = "_sites", 
                         ignore.chr = c("chrM", "MT"), ignore.string = "unique", ric = NULL,
                         expand.by.count = NULL, as.fragments = FALSE){

  if(!missing(path) & !missing(files)){
    stop("Both 'path' and 'files' cannot be defined.",
         call. = FALSE)
  }

  if(missing(path) & missing(files)){
    stop("Provide either a single path to the files of interest
         or a character vector of files to import.",
         call. = FALSE)
  }

  if(!missing(path)){
    ff <- list.files(path = path,
                     full.names = TRUE,
                     pattern = pattern)
  } else{
    ff <- files
  }

  if(!is.null(ric)) {
    if(dir.exists(ric)) {
      ric_files <- list.files(path = ric,
                              full.names = TRUE,
                              pattern = pattern)
      ff <- c(ff, ric_files)
    } else if(file.exists(ric)) {
      ff <- c(ff, ric)
    } else {
      warning("The specified 'ric' path does not exist. No random integration controls added.",
              call. = FALSE)
    }
  }

  ff <- ff[!grepl(ignore.string, ff)]

  sites <- lapply(
    X = ff,
    FUN = function(x){
      gr <- import(con = x,
                   format = sub(pattern=".", replacement = "", x = pattern))

      if(!is.null(levels.style)){
        seqlevelsStyle(gr) <- levels.style
        }

      gr <- gr[!seqnames(gr) %in% ignore.chr]

      if(!is.null(expand.by.count)) {
        if(expand.by.count %in% colnames(mcols(gr))) {
          counts <- round(mcols(gr)[[expand.by.count]])
          counts[is.na(counts) | counts < 1] <- 1
        } else {
          stop("Defined count column not found in sites.", call. = FALSE)
        }

        expanded.gr <- sapply(
          X = seq_along(gr),
          FUN = function(x){rep(gr[x], counts[x])}
        )

        gr <- unlist(as(expanded.gr, "GRangesList"))
        rm(expanded.gr)
      }

      if(!keep.metadata){
        mcols(gr) <- NULL
      }

      if(sort){
        gr <- sortSeqlevels(gr)
        gr <- sort(gr, ignore.strand = TRUE)
      }

      return(gr)
      }
    )

  # nn <- gsub( paste( ".*/(.*?)\\", pattern, sep = "" ), "\\1", ff )
  nn <- gsub(paste(".*/(.*?)\\", pattern, "(\\.\\w+)*(\\.\\w+)?$", sep = ""), "\\1", ff)

  if(!is.null(rm.string)){
    nn <- gsub(rm.string, "", nn)
  }

  if(any(duplicated(nn))){
    stop("Duplicated names found.",
         "\n",
         "The duplicated filenames are: ", nn[duplicated(nn)], ".",
         call. = FALSE)
  } else{
    names(sites) <- nn
  }

  if(!is.null(genome.obj)){
    levs <- lapply(X = sites,
                   FUN = function(x){
                     seqlevels(x)
                    }
    )

    levs <- unique(do.call(c, levs))

    out <- setdiff(levs, seqlevels(genome.obj))
    levs <- levs[levs %in% seqlevels(genome.obj)]

    if(length(out) > 0){
      warning("The following seqlevels were removed from integration site datasets: ",
              paste( out, collapse = ", "),
              call. = FALSE )
    }

    sites <- lapply(X = sites,
                    FUN = function(x){
                      seqlevels(x, pruning.mode = "coarse") <- levs
                      return(x)
                      }
                    )
    }

  if(as.fragments){
    return(sites)
  } else {
    sites <- SiteList(sites)
    return(sites)
  }
}