#' Estimate Logo Credible Intervals
#'
#' Estimates credible intervals for each base at each position in a sequence logo.
#'
#' @param sites A SiteList object.
#' @param genome.obj The BSgenome object for the genome of interest.
#' @param seq.len The desired length of the expanded sequences.
#' @param ignore.strand Boolean. Whether or not to ignore strandedness when extracting sequences.
#' Defaults to FALSE. If TRUE, the sequence from the forward strand is returned.
#' @param current.start The position in the target site duplication currently described by the start coordinates in sites.
#' This is used internally for centering the integration site coordinates.
#' @param tsd The total length of the target site duplication.
#' This is used along with current.start for centering the integration site coordinates.
#' @param alpha A scalar or vector of pseudo-counts to add to the observed counts for each nucleotide at each position.
#' Defaults to c(A = 0.3, C = 0.2, G = 0.2, T = 0.3). Be mindful of order!
#' @param ci A numeric vector denoting the upper- and lower-boundaries of the credibile interval. Defaults to c(0.025, 0.975).
#' @param seq.samp If not NULL, the number of sites to sample from a given dataset. This can significantly speed up processing time,
#' at the expense of including every site in the sequence logo. Defaults to NULL
#' @param dir.samp The number of samples to draw from the Dirichlet posterior distribution. Defaults to 5000.
#' @param seed The seed value for random sampling.
#' @param return.plot Boolean.
#' Whether or not to plot the results. Defaults to TRUE.
#' @param add.title Boolean. Whether or not to add the sample name to the plot. Defaults to FALSE.
#' @param wrap Boolean. Used with return.plot.
#' Whether or not to wrap the returned plots into a single output.
#' When TRUE, arguments ncol and nrow can be used to format the output.
#' @param ... Used to pass additional plotting parameters to patchwork::wrap_plots().
#' Examples include ncol and nrow.
#' 
#' @return If return.plot = TRUE, a list of ggplot2 objects or a single ggplot2 object containing sequence logos
#' for the mapped integration site datasets in sites.
#' If return.plot = FALSE, a list of data frames storing mean and median bit values, along with the estimated CIs.
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hs1)
#' data(sites2)
#' logo_ci(sites = sites2,
#'         genome.obj = Hsapiens)
#'
#' @import ggplot2
#' @import GenomicRanges
#' @import Biostrings
#' @import patchwork
#' @importFrom withr with_seed
#'
#' @export
#'
logo_ci <- function(sites, genome.obj, seq.len = 12, ignore.strand = FALSE, current.start = 1, 
                    tsd = 5, alpha = c(A = 0.3, C = 0.2, G = 0.2, T = 0.3), ci = c(0.025, 0.975), seq.samp = NULL, dir.samp = 5000, 
                    seed = NULL, return.plot = TRUE, add.title = FALSE, wrap = FALSE, ...) {
  
    if(!validObject(sites)){
    stop("sites is not a valid SiteListObject.",
         call. = FALSE)
  }

  if(seq.len %% 2 != 0){
    warning("seq.len cannot be odd. Subtracting 1 from seq.len.")
    seq.len <- seq.len - 1
    }
  
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed }, add = TRUE)
  } else {
    on.exit({ 
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(.Random.seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  expanded.sites <- expand_sites(
    sites = sites@sites,
    seq.len = seq.len,
    genome.obj = genome.obj,
    current.start = current.start,
    tsd = tsd
  )

  if(ignore.strand) {
    expanded.sites <- lapply(
      X = expanded.sites,
      FUN = function(x){
        strand(x) <- "+"
        return(x)
        })
  }

  if(!is.null(seq.samp)) {
    if(!is.null(seed)) {
      set.seed(seed)
    }

    expanded.sites <- lapply(
      X = expanded.sites,
      FUN = function(x){
        if (length(x) > n.samp){
          sample(x = x, size = n.samp, replace = FALSE)
        } else {x}})
  }

  seqs <- getSeq(
    x = genome.obj,
    names = as(expanded.sites, "GRangesList"),
    as.character = FALSE
    )

  mat.list <- lapply(
    X = seqs,
    FUN = function(x){
      consensusMatrix(x, baseOnly = TRUE)[seq(1,4),]
      }
    )

  center <- (tsd + 1) / 2

  if(tsd %% 2 == 1){
    mp <- (seq.len / 2)
    zp <- mp - (center - 1)
    offset <- 1
  } else{
    mp <- seq.len / 2
    zp <- mp - (floor(center) - 1)
    offset <- 0
  }

  p.labs <- seq(1, (seq.len - offset)) - zp
  
  cis <- lapply(
    X = seq_along(mat.list),
    FUN = function(x) {
      counts <- mat.list[[x]]
      # A, C, G, T
      n.bases <- nrow(counts)
      n.pos <- ncol(counts)
      
      if(length(alpha) == 1) {
        alpha <- rep(alpha, n.bases)
      }
      if(length(alpha) != n.bases){
        stop("alpha must be scalar or length n.bases")
      }
      
      IC <- matrix(NA_real_, nrow = dir.samp, ncol = n.pos)
      heights <- array(NA_real_, dim = c(dir.samp, n.pos, n.bases))
      
      log2n.bases <- log2(n.bases)
      
      if(!is.null(seed)) {
        set.seed(seed)
      }

      for(s in seq_len(dir.samp)) {
        gam <- matrix(
          rgamma(n = n.bases * n.pos, shape = as.vector(counts) + rep(alpha, times = n.pos)),
          nrow = n.bases, 
          ncol = n.pos
        )
        
        # Normalize each column to get Dirichlet samples
        p <- apply(
          X = gam, 
          MARGIN = 2, 
          FUN = function(x) x / sum(x)
        )
        
        # Shannon entropy (bits) for each position
        logp <- ifelse(p > 0, log2(p), 0)
        H <- -colSums(p * logp)

        Nseqs <- colSums(counts)
        en <- (1 / logb(2)) * (n.bases - 1)/(2 * Nseqs)  # Should match correction in ggseqlogo
        
        IC[s, ] <- pmax(log2n.bases - (H + en), 0)

        for (pos in 1:n.pos) {
          heights[s, pos, ] <- p[, pos] * IC[s, pos]
        }
      }
      
      # ic.median <- apply(IC, 2, median)
      # ic.mean   <- colMeans(IC)
      # ic.ci     <- apply(IC, 2, quantile, probs = ci)
      
      letter.median <- apply(heights, c(2, 3), median)
      letter.mean   <- apply(heights, c(2, 3), mean)
      letter.ci.lo  <- apply(heights, c(2, 3), quantile, probs = ci[1])
      letter.ci.hi  <- apply(heights, c(2, 3), quantile, probs = ci[2])

      dat <- data.frame(
        position = rep(p.labs, each = 4),  # Repeat each position 4 times for A, C, G, T
        letter = rep(c("A", "C", "G", "T"), times = n.pos),  # Letters repeated for each position
        bit.mean = as.vector(t(letter.mean)),  # Flatten the mean matrix
        bit.median = as.vector(t(letter.median)),  # Flatten the median matrix
        bit.ci.low = as.vector(t(letter.ci.lo)),  # Flatten the lower CI matrix
        bit.ci.high = as.vector(t(letter.ci.hi)),  # Flatten the upper CI matrix
        sample = rep(names(mat.list)[x], 4 * length(p.labs))
      )

      dat$letter <- factor(dat$letter, levels = c("A", "T", "C", "G"))
      dat$position <- as.character(dat$position)
      dat$position <- factor(dat$position, levels = unique(as.character(dat$position)))

      if(add.title) {
        plot.title <- names(mat.list)[x]
      }

      if(return.plot) {
        return(base_ci_plot(x = dat, add.title = add.title, title = plot.title))
      } else {
        return(dat)
      }
    }
  )

  if(!return.plot) {
    cis <- do.call(rbind, cis)
  } else {
    if(wrap) {
      cis <- wrap_plots(cis, ...)
    }
  }

  return(cis)
}

base_ci_plot <- function(x, add.title = FALSE, title = NULL) {
  pp <- ggplot(data = x) +
    aes(x = position, y = bit.mean, fill = letter) +
    geom_bar(stat = "identity", position = position_dodge(width=0.9), color = "black") + 
    geom_errorbar(aes(ymin = bit.ci.low, ymax = bit.ci.high), 
                  width = 0.2, 
                  color = "black",
                  position = position_dodge(width=0.9)) +
    scale_fill_manual(values = c('#109648', '#D62839', '#255C99', '#F7B32B')) +
    theme_bw() + 
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.key = element_blank(),
      legend.key.size = unit(0.5, "cm"),
      legend.position = "top",
      legend.box.margin = margin(-10,-10,-10,-10),
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    ) +
    labs(x = "Position", y = "Bits")

  if(add.title) {
    pp <- pp + labs(title = title)
  }

  return(pp)
}
