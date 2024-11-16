#' Simulate Random Data
#'
#' Generate random fragments from a defined genome under experimentally matched conditions.
#'
#'@param genome.obj The BSgenome object of interest.
#'@param re.sites A character vector of restriction enzyme recognition sequences.
#'@param cut.after An integer that defines the nucleotide after which the
#'DNA is cleaved in the restriction enzyme recognition sequence. Defaults to 1.
#'@param n.sites The number of sites to generate.
#'Can be either a single value or a vector of values.
#'@param mean The target mean fragment size. Only used for random fragmentation.
#'@param sd The target fragment size standard deviation.
#'Only used for random fragmentation.
#'@param min.width The minimum desired fragment/read length.
#'Defaults to 10 bp.
#'@param max.width The maximum acceptable fragment length (insert size).
#'Defaults to 1200 bp.
#'@param max.bp The maximum read length.
#'Defaults to 150 bp.
#'@param iterations The number of random datasets to generate.
#'If > 1, the number of values given for n.sites must be 1.
#'In this case, random datasets consisting of the same number of sites are generated for each iteration.
#'@param n.cores The number of cores to use.
#'@param write.ranges Boolean. Whether or not to write out the coordinates of the generated random fragments.
#'Defaults to FALSE.
#'@param write.reads Boolean. Whether or not to save the generated reads as R1 and R2 fasta files.
#'Defaults to TRUE.
#'@param prefix The file prefix appended to saved files.
#'@param directory.path The directory path where output files will be saved.
#'@param compress Boolean. Whether or not to compress the output FASTA files. Defaults to TRUE.
#'@param collapse Boolean. Whether or not to collapse the generated fragments for each iteration
#'into a single object. Defaults to TRUE.
#'@param to.chr.ends Boolean. Whether or not to treat chromosome ends as capable of generating
#'potentially mappable fragments. Defaults to TRUE.
#'@param paired Boolean. Whether or not to save paired-end fasta files.
#'@param U3 Boolean. Whether or not mapping is done from the U3 end (5' LTR relative to top-strand).
#'Defaults to FALSE.
#'
#'@return Default behavior will save paired R1 and R2 fasta files to disk that mimic NGS data.
#'The fragment coordinates used to generate those fasta files will be returned as a GRanges object.
#'
#'@examples
#'\donttest{
#'library(BSgenome.Hsapiens.UCSC.hs1)
#'set.seed(1)
#'
#'simulate_random(
#'  genome.obj = Hsapiens,
#'  re.sites = NULL,
#'  n.sites = 10,
#'  mean = 400,
#'  sd = 100,
#'  min.width = 14,
#'  max.distance = 900,
#'  max.bp = 150,
#'  iterations = 1,
#'  write.ranges = FALSE,
#'  write.reads = FALSE
#')
#'
#'simulate_random(
#'  genome.obj = Hsapiens,
#'  re.sites = "TTAA",
#'  cut.after = 2,
#'  n.sites = 10,
#'  min.width = 14,
#'  max.distance = 900,
#'  max.bp = 150,
#'  iterations = 1,
#'  write.ranges = FALSE,
#'  write.reads = FALSE
#')
#'}
#'
#'@import methods
#'@import Biostrings
#'@import parallel
#'
#'@export
#'
simulate_random <- function( genome.obj,
                             re.sites = NULL,
                             cut.after = 1,
                             n.sites,
                             mean = 400,
                             sd = 100,
                             min.width = 10,
                             max.width = 1200,
                             max.bp = 150,
                             iterations = 1,
                             n.cores = 1,
                             write.ranges = FALSE,
                             write.reads = TRUE,
                             prefix = NULL,
                             directory.path = ".",
                             compress = TRUE,
                             paired = TRUE,
                             collapse = TRUE,
                             to.chr.ends = TRUE,
                             U3 = FALSE ){

  if( iterations > 1 & length( n.sites ) > 1 ){
    stop( "iterations and length(n.sites) cannot both be > 1.",
          call. = FALSE )
  }

  if( iterations > 1 & length( n.sites ) == 1 ){
    n.sites <- rep( n.sites, iterations )
  }

  message( "\n" )
  message( "Getting chromosome sequences...", "\n" )

  chr.seqs <- get_chromosome_seqs( genome.obj = genome.obj )

  if ( is.null( re.sites ) ){
    random <- TRUE
    re.cuts <- NULL
  } else{
    message( "Finding restriction enzyme target sequences...", "\n" )
    re.cuts <- digest( seqs.list = chr.seqs,
                       re.sites = re.sites,
                       cut.after = cut.after )
    random <- FALSE
  }

  if ( isTRUE( random ) ){
    if ( !is.numeric( mean ) | !is.numeric( sd ) )
      stop( "The mean and sd of the sampling distribution must be
            numeric arguments.",
            call. = FALSE )
  }

  frags <- mclapply( X = seq_along( n.sites ),
                     FUN = function(x){

                       message( "Set ", x )
                       message( "-----", "\n" )

                       tmp.n <- n.sites[x]

                       message( "Generating random site positions...", "\n" )

                       rand.sites <- random_sites( n.sites = tmp.n,
                                                   genome.obj = genome.obj )

                       message( "Getting fragment ranges...", "\n" )

                       rand.fragments <- make_fragments( insert.sites = rand.sites,
                                                         frag.sites = re.cuts,
                                                         random = random,
                                                         genome.obj = genome.obj,
                                                         mean = mean,
                                                         sd = sd,
                                                         to.chr.ends = to.chr.ends,
                                                         U3 = U3 )

                       if( write.reads ){
                         message( "Extracting fragment sequences...", "\n" )

                         frag.seqs <- getSeq( x = genome.obj,
                                              names = rand.fragments,
                                              as.character = FALSE )

                         message( "Trimming fragment sequences...", "\n" )

                         frag.trim <- trim_seqs( fragments = frag.seqs,
                                                 min.width = min.width,
                                                 max.width = max.width,
                                                 max.bp = max.bp )

                         message( "Saving FASTA files...", "\n" )

                         if ( is.null( directory.path ) ){
                           stop( "directory.path and prefix must be defined to save FASTA files.",
                                 call. = FALSE)
                         }

                         if( is.null( prefix ) ){
                           prefix <- paste0( "set_", x )
                         } else{
                           prefix <- paste0( prefix, "_", x )
                         }

                         save_fasta( reads = frag.trim,
                                     directory.path = directory.path,
                                     prefix = prefix,
                                     paired = paired,
                                     compress = compress )
                       }

                       return( rand.fragments )

                       },

                     mc.cores = n.cores

                     )

  if( is.null( prefix ) ){
    prefix <- "set"
  }

  if( isTRUE( collapse ) ){

    generated.fragments <- do.call( c, frags )

  } else{

    generated.fragments <- frags
    names( generated.fragments ) <- paste0( prefix, "_", seq_along( n.sites ) )

  }

  if( isTRUE( write.ranges ) ){

    message( "Writing generated fragment ranges...", "\n" )

    save( generated.fragments,
          file = paste0( directory.path, "/", prefix, "_generated_fragments.RData.gz" ),
          compression_level = 6 )

    message( "Done!", "\n" )

    } else {

      message( "Done!", "\n" )

      return( generated.fragments )

    }
  }
