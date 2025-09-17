#' Simulate Random Data
#'
#' Generate random fragments from a defined genome under experimentally matched conditions.
#'
#' @param genome.obj The BSgenome object of interest.
#' @param re.sites A character vector of restriction enzyme recognition sequences.
#' @param cut.after An integer that defines the nucleotide after which the
#' DNA is cleaved in the restriction enzyme recognition sequence. Defaults to 1.
#' @param n.sites The number of sites to generate.
#' Can be either a single value or a vector of values.
#' @param mean The target mean fragment size. Only used for random fragmentation.
#' @param sd The target fragment size standard deviation.
#' Only used for random fragmentation.
#' @param min.width The minimum desired fragment/read length.
#' Defaults to 25 bp.
#' @param max.width The maximum acceptable fragment length (insert size).
#' Defaults to 1200 bp.
#' @param max.bp The maximum read length.
#' Defaults to 150 bp.
#' @param iterations The number of random datasets to generate.
#' If > 1, the number of values given for n.sites must be 1.
#' In this case, random datasets consisting of the same number of sites are generated for each iteration.
#' @param n.cores The number of cores to use.
#' @param write.ranges Boolean. Whether or not to write out the coordinates of the generated random fragments.
#' Defaults to FALSE.
#' @param write.reads Boolean. Whether or not to save the generated reads as R1 and R2 fastq files.
#' Defaults to TRUE.
#' @param prefix The file prefix appended to saved files. Can be NULL, a single character value, or a vector of character values.
#' @param directory.path The directory path where output files will be saved.
#' @param compress Boolean. Whether or not to compress the output FASTA files. Defaults to TRUE.
#' @param collapse Boolean. Whether or not to collapse the generated fragments for each iteration
#' into a single object. Defaults to TRUE.
#' @param to.chr.ends Boolean. Whether or not to treat chromosome ends as capable of generating
#' potentially mappable fragments. Defaults to TRUE.
#' @param paired Boolean. Whether or not to save paired-end fasta files. Defaults to TRUE.
#' @param U3 Boolean. Whether or not mapping is done from the U3 end (5' LTR relative to top-strand).
#' Defaults to FALSE.
#' @param ltr.seq The LTR sequence sequenced in the experiment, given in the 5' to 3' orientation.
#' Should be a character string or NULL. Must be provided.
#' @param linker.seq The linker sequence sequenced in the experiment, given in the 5' to 3' orientation.
#' Should be a character string or NULL. Must be provided.
#' @param mismatch.prob The NGS mismatch rate. Defaults to 0.0025.
#' @param indel.prob The NGS indel rate. Defaults to 2.5E-5.
#' @param base.score The quality score value to assign to all bases. Offset is added to this value. Default 40.
#' @param offset The ASCII offset of the quality scores. Defaults to 33.
#' @param n.clones The number of clonal sites to generate. Defaults to 0. Can be given as a list of vectors equal in length to
#' the prefix vector or the number of iterations.
#' @param clone.fracs A numeric vector defining the fraction of sites occupied by each clone. Defaults to 0. Can be given as a list of
#' vectors equal in length to the prefix vector or the number of iterations.
#' Must be of length == n.clones, all values must be between 0 and 1, and the sum of values must be <= 1.
#' @param as.fasta Boolean. Whether or not to save the generated reads in FASTA format. If FALSE, reads are saved as FASTQ.
#' Defaults to FALSE.
#' @param verbose Boolean. Whether or not to print status updates. Defaults to FALSE.
#' @param chr.omit A character vector of chromosomes to omit. Defaults to c("chrM", "MT")
#' @param noisy.scores Boolean. Whether or not to introduce simulated noise into quality scores. Defaults to FALSE.
#' @param read.name.prefix The naming prefix for read names. Defaults to "fragment_".
#' @param clone.chroms A vector denoting the chromosome harboring each clonal site. Defaults to NULL.
#' If not NULL, must be of length equal to n.clones.
#' @param clone.positions A vector denoting the position of each clonal site. Defaults to NULL.
#' If not NULL, must be of length equal to n.clones.
#' @param clone.strands A vector denoting the strand of each clonal site. Defaults to NULL.
#' If not NULL, must be of length equal to n.clones.
#' @param use.frag.names Boolean. Whether or not use names from the input fragment GRanges for sequence names. Defaults to TRUE.
#'
#' @return Default behavior will save paired R1 and R2 fastq files to disk.
#' The fragment coordinates used to generate those fastq files will be returned as a GRanges object.
#'
#' @examples
#' \donttest{
#' library(BSgenome.Hsapiens.UCSC.hs1)
#' set.seed(1)
#'
#' simulate_random(
#'   genome.obj = Hsapiens,
#'   re.sites = NULL,
#'   n.sites = 10,
#'   mean = 400,
#'   sd = 100,
#'   min.width = 25,
#'   max.distance = 1200,
#'   max.bp = 150,
#'   iterations = 1,
#'   write.ranges = FALSE,
#'   write.reads = FALSE
#' )
#'
#' simulate_random(
#'   genome.obj = Hsapiens,
#'   re.sites = "TTAA",
#'   cut.after = 2,
#'   n.sites = 10,
#'   min.width = 25,
#'   max.distance = 1200,
#'   max.bp = 150,
#'   iterations = 1,
#'   write.ranges = FALSE,
#'   write.reads = FALSE
#' )
#'}
#'
#' @import methods
#' @import Biostrings
#' @import parallel
#'
#' @export
#'
simulate_random <- function(genome.obj,
                            re.sites = NULL,
                            cut.after = 1,
                            n.sites,
                            mean = 400,
                            sd = 100,
                            min.width = 25,
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
                            collapse = FALSE,
                            to.chr.ends = TRUE,
                            U3 = FALSE,
                            ltr.seq,
                            linker.seq,
                            mismatch.prob = 0.0025,
                            indel.prob = 2.5E-5,
                            offset = 33,
                            as.fasta = FALSE,
                            n.clones = 0,
                            clone.fracs = 0,
                            verbose = FALSE,
                            chr.omit = c("chrM", "MT"),
                            base.score = 40,
                            noisy.scores = FALSE,
                            read.name.prefix = "fragment_",
                            use.frag.names = TRUE,
                            clone.chroms = NULL,
                            clone.positions = NULL,
                            clone.strands = NULL){

  #TO-DO: replace mclapply with parLapply or BiocParallel
  #TO-DO: Make sure fragments aren't generated for chr.omit chromosomes

  print_message <- function(...) {
    if(isTRUE(verbose)) {
      message(...)
      }
  }

  if(iterations > 1 && length(n.sites) > 1){
    stop("iterations and length(n.sites) cannot both be > 1.",
          call. = FALSE)
  }

  if(iterations > 1 && length(n.sites) == 1){
    n.sites <- rep(n.sites, iterations)
  }

  if(length(prefix) > 1 && length(n.sites) > 1){
    if(length(prefix) != length(n.sites)){
      stop("When length(prefix) > 1, length(prefix) must equal length(n.sites) or iterations.",
           call. = FALSE)
    }
  }

  if(length(n.sites) == 1 && length(prefix) > 1){
    n.sites <- rep(n.sites, length(prefix))
  }

  if(all(is.list(n.clones), is.list(clone.fracs))){
    if(length(n.clones) != length(clone.fracs)){
      stop("length(n.clones) must equal length(clone.fracs)")
    }
    if((length(n.clones) != length(n.sites))){
      stop("length(n.clones) must be the same as the number of datasets being generated
           (defined by n.sites, iterations, and/or prefix).",
           call. = FALSE)
    }
  } else if(!all(!is.list(n.clones), !is.list(clone.fracs))){
    stop("n.clones and clone.fracs must both be lists or not be lists.",
         call. = FALSE)
  }

  print_message("\n")
  print_message("Getting chromosome sequences...", "\n")

  if(is.null(re.sites)){
    random <- TRUE
    re.cuts <- NULL
  } else{
    print_message("Finding restriction enzyme target sequences...", "\n")
    chr.seqs <- setNames(
      lapply(
        X = names(genome.obj),
        FUN = function(chr){getSeq(genome.obj, chr)}
      ),
      names(genome.obj)
    )

    chr.seqs <- chr.seqs[!names(chr.seqs) %in% chr.omit]

    re.cuts <- digest(
      chr.seqs = chr.seqs,
      re.sites = re.sites,
      cut.after = cut.after
      )
    random <- FALSE
  }

  if(isTRUE(random)){
    if(!is.numeric(mean) | !is.numeric(sd))
      stop("The mean and sd of the sampling distribution must be numeric.",
            call. = FALSE)
    }

  frags <- mclapply(
    X = seq_along(n.sites),
    FUN = function(x){
      print_message("Set ", x)
      print_message("-----", "\n")

      tmp.n <- n.sites[x]

      if(is.list(n.clones)){
        tmp.n.clones <- n.clones[[x]]
      } else{
        tmp.n.clones <- n.clones
      }

      if(is.list(clone.fracs)){
        tmp.clone.fracs <- clone.fracs[[x]]
      } else{
        tmp.clone.fracs <- clone.fracs
      }

      print_message("Generating random site positions...", "\n")

      rand.sites <- random_sites(
        n.sites = tmp.n,
        genome.obj = genome.obj,
        n.clones = tmp.n.clones,
        clone.fracs = tmp.clone.fracs,
        chr.omit = chr.omit,
        clone.chroms = clone.chroms,
        clone.positions = clone.positions,
        clone.strands = clone.strands
        )

      print_message("Getting fragment ranges...", "\n")

      rand.fragments <- make_fragments(
        insert.sites = rand.sites,
        frag.sites = re.cuts,
        random = random,
        genome.obj = genome.obj,
        mean = mean,
        sd = sd,
        to.chr.ends = to.chr.ends,
        U3 = U3,
        read.name.prefix = read.name.prefix
        )

      if(write.reads){
        print_message("Extracting fragment sequences...", "\n")

        batch.size <- min(tmp.n, 20000)
        if(batch.size < 20000){
          n.batches <- 1
        } else{
          n.batches <- ceiling(tmp.n / batch.size)
        }

        name.idx.start <- 0
        print_message("Trimming fragment sequences...", "\n")

        batch.starts <- cumsum(c(0, rep(batch.size, n.batches - 1)))
        
        batch.results <- mclapply(
          X = seq_len(n.batches),
          FUN = function(batch) {
            start.idx <- (batch - 1) * batch.size + 1
            end.idx <- min(batch * batch.size, tmp.n)
            batch.frags <- rand.fragments[start.idx:end.idx]

            batch.seqs <- getSeq(
              x = genome.obj,
              names = batch.frags,
              as.character = FALSE
            )

            if (use.frag.names) {
              names(batch.seqs) <- batch.frags$name
            }

            batch.trim <- trim_seqs(
              fragments = batch.seqs,
              min.width = min.width,
              max.width = max.width,
              max.bp = max.bp,
              ltr.seq = ltr.seq,
              linker.seq = linker.seq,
              mismatch.prob = mismatch.prob,
              indel.prob = indel.prob,
              offset = offset,
              base.score = base.score,
              as.fasta = as.fasta,
              noisy.scores = noisy.scores,
              read.name.prefix = read.name.prefix,
              name.idx.start = batch.starts[batch],
              use.frag.names = use.frag.names,
              paired = paired
            )

            list(batch.trim = batch.trim, batch.size = length(batch.trim$R1))
          },
          mc.cores = ifelse((length(n.sites) == 1 & iterations == 1), n.cores, 1)
        )

        if(as.fasta){
          frag.trim <- list(
            R1 = do.call(c, lapply(batch.results, function(res) res$batch.trim$R1))
          )
          if(paired) {
            frag.trim$R2 <- do.call(c, lapply(batch.results, function(res) res$batch.trim$R2))
          }
        } else {
          frag.trim <- list(
            R1 = list(
              seq = unlist(lapply(batch.results, function(res) res$batch.trim$R1$seq)),
              qual = unlist(lapply(batch.results, function(res) res$batch.trim$R1$qual))
            )
          )
          if(paired) {
            frag.trim$R2 <- list(
              seq = unlist(lapply(batch.results, function(res) res$batch.trim$R2$seq)),
              qual = unlist(lapply(batch.results, function(res) res$batch.trim$R2$qual))
            )
          }
        }

        print_message("Saving output files...", "\n")

        if (is.null(directory.path)){
          stop("directory.path and prefix must be defined to save files.",
               call. = FALSE)
          }

        if(is.null(prefix)){
          tmp.prefix <- paste0("set_", x)
        } else if(length(prefix) > 1){
          tmp.prefix <- prefix[x]
        } else{
          tmp.prefix <- paste0(prefix, "_", x)
        }

        if(isTRUE(as.fasta)){
          save_fasta(
            reads = frag.trim,
            directory.path = directory.path,
            prefix = tmp.prefix,
            paired = paired,
            compress = compress
          )
        } else{
          save_fastq(
            reads = frag.trim,
            directory.path = directory.path,
            prefix = tmp.prefix,
            paired = paired,
            compress = compress
            )
          }
        }
      return(rand.fragments)
      },
    mc.cores = ifelse((length(n.sites) > 1 & iterations > 1), n.cores, 1)
    )

  if(is.null(prefix)){
    prefix <- "set"
  }

  if(isTRUE(collapse)){
    generated.fragments <- do.call(c, frags)
  } else if(length(prefix) > 1){
    generated.fragments <- frags
    names(generated.fragments) <- prefix[seq_along(prefix)]
  } else{
    generated.fragments <- frags
    names(generated.fragments) <- paste0(prefix, "_", seq_along(n.sites))
  }

  if(isTRUE(write.ranges)){
    print_message("Writing generated fragment ranges...", "\n")

    if(length(prefix) == 1 && !isTRUE(collapse)){
      lapply(
        X = seq_along(generated.fragments),
        FUN = function(x){
          assign(paste0(prefix, "_", x), generated.fragments[[x]])
          save(
            list = paste0(prefix, "_", x),
            file = paste0(directory.path, "/", prefix, "_", x, "_generated_fragments.RData.gz"),
            compression_level = 6
          )
        }
      )
    } else if(isTRUE(collapse)){
      save(
        generated.fragments,
        file = paste0(directory.path, "/", "generated_fragments.RData.gz"),
        compression_level = 6
      )
    } else if(length(prefix) > 1 && !isTRUE(collapse)){
      lapply(
        X = seq_along(prefix),
        FUN = function(x){
          tmp.name <- paste0(prefix[x])
          assign(tmp.name, generated.fragments[[x]])
          save(
            list = tmp.name,
            file = paste0(directory.path, "/", prefix[x], "_generated_fragments.RData.gz"),
            compression_level = 6
          )
        }
      )
    }
    print_message("Done!", "\n")
    } else {
      print_message("Done!", "\n")
      return(generated.fragments)
    }
  }
