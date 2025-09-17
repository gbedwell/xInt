#' Plot by Chromosome
#'
#' Generates karyotype plots for genomic regions and features across chromosomes.
#'
#' @param sites A SiteList object.
#' @param feature.list A named list of genomic features to plotted.
#' @param genome.obj A BSgenome object.
#' @param conditions Optional. A vector of conditions.
#' Must be the same length as sites. 
#' When given, the datasets belonging to each condition are collapsed into a pooled dataset.
#' @param win.size The window size for feature density calculations. Defaults to 1E6.
#' @param stack.cutoff Defines the maximum number of layers each plotting level can contain.
#' Defaults to 4.
#' @param chr.omit A character vector of chromosomes to omit from the plot. Defaults to c("chrM", "MT").
#' Must be NULL if chr.include is given.
#' @param chr.include A character vector of chromosomes to include in the plot. Defaults to NULL.
#' @param cond.omit A vector of conditions to omit from the plot.
#' @param label.size The size of the feature labels. Defaults to 0.8.
#'
#' @import karyoploteR
#'
#' @return A list of plots, invisibly returned.
#'
#' @examples
#' plot_by_chrom(sites, feature.list = list(genes = genes), genome.obj = Hsapiens)
#'
#' @export
plot_by_chrom <- function(sites, feature.list, genome.obj, conditions = NULL,
                          win.size = 1E6, stack.cutoff = 4, chr.omit = NULL,
                          chr.include = NULL, cond.omit = NULL, label.size = 0.8) {
  
  if (!validObject(sites)) {
    stop("sites is not a valid SiteList object.", call. = FALSE)
  }
  
  if (!is.list(feature.list) || is.null(names(feature.list))) {
    stop("feature.list must be a named list.",
         call. = FALSE)
  }
  
  if (!is.null(conditions) && length(conditions) != length(sites)) {
    stop("conditions must have the same length as sites.", call. = FALSE)
  }
  
  if (!is.null(conditions)) {
    cond.list <- lapply(
      X = unique(conditions),
      FUN = function(x) {
        which(conditions == x)
      }
    )
    
    names(cond.list) <- unique(conditions)
    
    sites <- collapse_sites(sites, cond.list, sorted = TRUE)
    conditions <- unique(conditions)
  } else {
    conditions <- names(sites)
  }
  
  if (!is.null(cond.omit)) {
    sites <- sites[!names(sites) %in% cond.omit]
  }
  
  if(!is.null(chr.omit) && !is.null(chr.include)) {
    stop("Values given for both chr.omit and chr.include.", 
         call. = FALSE)
  }
  
  if(!is.null(chr.omit)) {
    target.chroms <- seqnames(genome.obj)[!seqnames(genome.obj) %in% chr.omit]
  } else if(!is.null(chr.include)) {
    target.chroms <- seqnames(genome.obj)[seqnames(genome.obj) %in% chr.include]
  } else {
    target.chroms <- seqnames(genome.obj)
  }
  
  feat.colors <- get_custom_palette(length(feature.list) + 1)
  
  region.lists <- lapply(
    X = seq_along(sites),
    FUN = function(site.idx) {
      combo.list <- list(
        list(
          data = sites[[1]], 
          color = feat.colors[1],
          name = "ISs"
          )
      )
      
      for(feat.idx in seq_len(length(feat.colors) - 1)) {
        combo.list[[feat.idx + 1]] <- list(
          data = feature.list[[feat.idx]],
          color = feat.colors[feat.idx + 1],
          name = names(feature.list)[feat.idx]
        )
      }
      
      return(combo.list)
    }
  )
  
  kp <- plotKaryotype(
    genome = Hsapiens, 
    chromosomes = target.chroms,
    ideogram.plotter = NULL,
    cytobands = GRanges()
  )
  
  r0 <- 0
  r1.increment <- 0.125
  
  plots <- lapply(
    X = region.lists,
    FUN = function(ds) {
      for(feat.idx in seq_along(ds)) {
        if(feat.idx == 1) {
          kpPlotRegions(
            kp, 
            data = ds[[feat.idx]]$data, 
            col = transparent(ds[[feat.idx]]$color, amount = 0.5), 
            data.panel = "ideogram"
          )
          
          kpPlotDensity(
            kp, 
            data = ds[[feat.idx]]$data, 
            data.panel = "ideogram", 
            col = transparent(ds[[feat.idx]]$color, amount = 0.5), 
            window.size = win.size
          )
          kpAddLabels(
            kp, 
            labels = ds[[feat.idx]]$name, 
            label.margin = 0.01,  
            side = "right", 
            data.panel = "ideogram",
            cex = label.size
          )
        } else {
          bins <- disjointBins(ds[[feat.idx]]$data)
          kpPlotRegions(
            kp, 
            data = ds[[feat.idx]]$data, 
            num.layers = ifelse(max(bins) > stack.cutoff, stack.cutoff, 2),
            col = transparent(ds[[feat.idx]]$color, amount = 0.5), 
            r0 = r0, 
            r1 = r0 + r1.increment
          )
          kpPlotDensity(
            kp, 
            data = ds[[feat.idx]]$data, 
            data.panel = 1, 
            col = transparent(ds[[feat.idx]]$color, amount = 0.7), 
            window.size = win.size, 
            r0 = r0, 
            r1 = r0 + 2 * r1.increment
          )
          kpAddLabels(
            kp, 
            labels = ds[[feat.idx]]$name, 
            label.margin = 0.01,  
            side = "right", 
            data.panel = 1,
            r0 = r0, 
            r1 = r0 + 2 * r1.increment,
            cex = label.size
          )
          
          r0 <- r0 + 2 * r1.increment
        }
      }
    }
  )
}
