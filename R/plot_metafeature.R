#' Plot Meta-Feature Distributions
#'
#' Plots meta-feature distributions using the output of metafeature().
#'
#' @param x The output of metafeature().
#' @param bins The number of bins to use for binning the data. Defaults to 100.
#' @param order The desired order of the data in the figure legend. All items must be present in x$conditions.
#' @param average.samples Boolean. Whether or not to plot the mean and standard deviation across replicates.
#' Defaults to TRUE. When FALSE, replicate data is pooled.
#' @param shade.under A character vector of length 1 defining a reference condition to shade under.
#' Usually used with a random integration control.
#' Helps highlight deviations from the reference.
#' @param bp.range A numeric vector of length 2 defining the basepair range to display (e.g., c(-1000, 1000)).
#' If NULL, uses relative position (0 to 1).
#' By setting this, you are implying that all features are the same length.
#' @param y.scale One of "percent" (default) or "fraction". Determines the scale of the y-axis.
#' 
#' @return A ggplot2 object for plotting meta-feature distributions.
#'
#' @import ggplot2
#'
#' @export
#'
plot_metafeature <- function(x, bins = 100, order = NULL, average.samples = TRUE, shade.under = NULL, 
                             bp.range = NULL, y.scale = c("percent", "fraction")) {
  if(!all(c("rel.position", "dataset") %in% colnames(x))) {
    stop("Input does not look like the output of metafeature().", call. = FALSE)
  }
    
  if(!is.null(order) && !all(order %in% unique(x$condition))) {
    stop("Not all items in order are present in the input data.", call. = FALSE)
  }

  y.scale = match.arg(y.scale)
  
  split.df <- split(x, f = x$condition)
  
  bin_data <- function(values, n.bins) {
    breaks <- seq(0, 1, length.out = n.bins + 1)
    bin.centers <- (breaks[-1] + breaks[-length(breaks)]) / 2
    bins <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    # Return the bin centers for each value
    return(bin.centers[bins])
  }
  
  bin.df <- lapply(
    X = split.df,
    FUN = function(dat) {
      if(average.samples && length(unique(dat$dataset)) > 1) {
        # Split by dataset within this condition
        dataset.splits <- split(dat, f = dat$dataset)
        
        # Process each dataset separately
        dataset.bins <- lapply(dataset.splits, function(ds) {
          bin.pos <- bin_data(ds$rel.position, n.bins = bins)
          bin.counts <- table(bin.pos)
          bin.pos <- as.numeric(names(bin.counts))
          bin.counts <- as.numeric(bin.counts)

          if(y.scale == "percent") {
            bin.frac <- bin.counts / sum(bin.counts) * 100
          } else {
            bin.frac <- bin.counts / sum(bin.counts)
          }
          
          data.frame(pos = bin.pos, frac = bin.frac)
        })
        
        # Combine and average across datasets
        all.pos <- sort(unique(unlist(lapply(dataset.bins, function(db) db$pos))))
        
        # Calculate mean and standard deviation for each position
        avg.fracs <- sapply(all.pos, function(p) {
          fracs <- sapply(dataset.bins, function(db) {
            idx <- which(db$pos == p)
            if(length(idx) > 0) db$frac[idx] else 0
          })
          mean(fracs)
        })
        
        sd.fracs <- sapply(all.pos, function(p) {
          fracs <- sapply(dataset.bins, function(db) {
            idx <- which(db$pos == p)
            if(length(idx) > 0) db$frac[idx] else 0
          })
          if(length(fracs) > 1) sd(fracs) else 0
        })
        
        bin.dat <- data.frame(
          pos = all.pos,
          frac = avg.fracs,
          frac.sd = sd.fracs,
          id = dat$condition[1]
        )
      } else {
        bin.pos <- bin_data(dat$rel.position, n.bins = bins)
        bin.counts <- table(bin.pos)
        bin.pos <- as.numeric(names(bin.counts))
        bin.counts <- as.numeric(bin.counts)
        bin.frac <- bin.counts / sum(bin.counts)
        bin.dat <- data.frame(
          pos = bin.pos,
          frac = bin.frac,
          frac.sd = 0,
          id = dat$condition[1]
        )
      }
      return(bin.dat)
    }
  )
  
  bin.df <- do.call(rbind, bin.df)
  rownames(bin.df) <- NULL

  if(!is.null(bp.range)) {
    min.bp <- bp.range[1]
    max.bp <- bp.range[2]
    bin.df$pos <- min.bp + bin.df$pos * (max.bp - min.bp)
    x.label <- "Distance (bp)"
  } else {
    x.label <- "Relative Position"
  }
  
  unique.ids <- unique(bin.df$id)
  if (!is.null(shade.under) && shade.under %in% unique.ids) {
    other.ids <- unique.ids[unique.ids != shade.under]
    if (!is.null(order)) {
      other.ids <- intersect(order, other.ids)
    }
    ordered.ids <- c(other.ids, shade.under)
    colors <- c(get_custom_palette(length(other.ids)), "black")
    names(colors) <- ordered.ids
    
    # Reorder the data
    bin.df$id <- factor(bin.df$id, levels = ordered.ids)
  } else {
    if(!is.null(order)) {
      ordered.ids <- intersect(order, unique.ids)
    } else {
      ordered.ids <- unique.ids
    }
    
    colors <- get_custom_palette(length(ordered.ids))
    names(colors) <- ordered.ids
    
    bin.df$id <- factor(bin.df$id, levels = ordered.ids)
  }
  
  bin.plot <- ggplot(bin.df) +
    aes(x = pos, y = frac, color = id)
  
  if(!is.null(shade.under)) {
    ref.data <- bin.df[bin.df$id == shade.under, ]
    bin.plot <- bin.plot +
      geom_ribbon(data = ref.data, 
                  aes(x = pos, ymin = 0, ymax = frac),
                  alpha = 0.3, color = NA, fill = "gray50", inherit.aes = FALSE) +
      scale_fill_discrete(guide = "none")
  }
  
  if(y.scale == "percent") {
    y.label <- "Percent Integration"
  } else {
    y.label <- "Fraction Integration"
  }

  bin.plot <- bin.plot +
    geom_line(linewidth = 1) +
    geom_ribbon(
      aes(y = frac, ymin = frac - frac.sd, ymax = frac + frac.sd, fill = id),
      alpha = .2,
      linewidth = 0
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_bw() + 
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.key = element_blank(),
      legend.position = "top",
      legend.box.margin = margin(-10,-10,-10,-10),
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    ) +
    guides(
      color = guide_legend(override.aes = list(fill = NA)),
      fill = "none"
    ) + 
    # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    # scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    labs(
      x = x.label,
      y = y.label
    )
  
  return(bin.plot)
}