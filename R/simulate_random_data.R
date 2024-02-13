#' Simulate Random Data
#'
#' Wrapper function that automatically calls the functions used to generate random genomic fragments. If re.sites is NULL, random fragmentation is assumed. When random fragmentation is assumed, mean and sd must be defined in the function call (default values in make_fragments() are not used). When re.sites is not NULL, the function tries to find re.sites in the chromosome sequences and uses these positions to make fragments. See the individual functions for more information. This function does not use any of the parallelized functions.
#'
#'@param genome.obj The name of the genome object of interest.
#'@param re.sites A character vector of restriction enzyme recognition sequences.
#'@param cut.after A single integer that defines the nucleotide after which the DNA is cleaved. Defaults to 1.
#'@param n.sites The number of sites to generate.
#'@param mean The target mean fragment size. Only used for random fragmentation.
#'@param sd The target fragment size standard deviation. Only used for random fragmentation.
#'@param min.width The minimum desired fragment/read length. Default 14 bp.
#'@param max.distance The maximum desired inner distance. Default 1000 bp.
#'@param max.bp The maximum read length. Default 150 bp.
#'@param iterations The number of random fragment datasets to generate.
#'@param n.cores The number of cores to use. Enables parallelization.
#'@param write.ranges Logical. Whether or not to write out the coordinates of the generated random fragments. Defaults to FALSE.
#'@param directory.path The parent directory path.
#'@param compress Boolean. Whether or not to compress the output FASTA files. Defaults to TRUE.
#'@param collapse Boolean. Whether or not to collapse the generated fragments for each iteration into a single object. Defaults to TRUE.
#'
#'@import BSgenome
#'@import parallel
#'
#'@export
#'
simulate_random_data <- function( genome.obj,
                                  re.sites = NULL,
                                  cut.after = 1,
                                  n.sites,
                                  mean,
                                  sd,
                                  min.width = 14,
                                  max.distance = 1000,
                                  max.bp = 150,
                                  iterations = 1,
                                  n.cores = 1,
                                  write.ranges = FALSE,
                                  prefix = NULL,
                                  directory.path = NULL,
                                  compress = TRUE,
                                  collapse = TRUE ){

  cat( "\n" )
  cat( "Getting chromosome sequences...", "\n\n" )

  chr.seqs <- get_chromosome_seqs( genome.obj = genome.obj )

  if ( is.null( re.sites ) ){
    random <- TRUE
    re.cuts <- NULL
  } else{
    cat( "Finding restriction enzyme target sequences...", "\n\n" )
    re.cuts <- digest( string.list = chr.seqs,
                       re.sites = re.sites,
                       cut.after = cut.after )
    random <- FALSE
  }

  if ( isTRUE( random ) ){
    if ( missing( mean ) | missing( sd ) )
      stop( "The mean and sd of the sampling distribution must be defined to generate fragments via random fragmentation.",
            call. = FALSE )
  }

  frags <- mclapply( X = 1:iterations,
                     FUN = function(x){

                       cat( "Set", x, "\n" )
                       cat( "-----", "\n\n" )

                       cat( "Generating random site positions...", "\n\n" )

                       rand.sites <- random_sites( n.sites = n.sites,
                                                   genome.obj = genome.obj )

                       cat( "Getting fragment ranges...", "\n\n" )

                       rand.fragments <- make_fragments( insert.sites = rand.sites,
                                                         frag.sites = re.cuts,
                                                         random = random,
                                                         genome.obj = genome.obj,
                                                         mean = mean,
                                                         sd = sd )

                       cat( "Extracting fragment sequences...", "\n\n" )

                       frag.seqs <- getSeq( x = genome.obj,
                                            names = rand.fragments,
                                            as.character = FALSE )

                       cat( "Trimming fragment sequences...", "\n\n" )

                       frag.trim <- trim_seqs( fragments = frag.seqs,
                                               min.width = min.width,
                                               max.distance = max.distance,
                                               max.bp = max.bp )

                       cat( "Saving FASTA files...", "\n\n" )

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
                                   compress = compress )

                       return( rand.fragments )

                       },

                     mc.cores = n.cores

                     )

  if( isTRUE( collapse ) ){

    generated.fragments <- do.call( c, frags )

  } else{

    generated.fragments <- frags

  }

  if( isTRUE( write.ranges ) ){

    cat( "Writing generated fragment ranges...", "\n\n" )

    save( generated.fragments,
          file = paste0( directory.path, "/", "generated_fragments.RData.gz" ),
          compression_level = 6 )

    cat( "Done!", "\n\n" )

    } else {

      cat( "Done!", "\n\n" )

      return( generated.fragments )

    }

  }
