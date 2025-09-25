#' Plot Feature-Set Overlap
#'
#' Plots global feature-set overlap for each condition in an xIntOverlap object. 
#'
#' @param xint.obj An xIntOverlap object.
#' @param stats.df A data frame containing significance values for each comparison.
#' Must have the columns 'comparison' and 'p.adj' (even if the p-values are not adjusted).
#' The outputs of e.g., global_comparisons(), overlap_comparisons() with type = 'global',
#' and overlap_mean_tests() will all work.
#' @param numeric.p Boolean. When FALSE (default), significance is denoted symbolically.
#' p > 0.05 is 'ns', 0.01 >= p < 0.05 is '*', 0.001 >= p < 0.01 is '**', and p < 0.001 is '***'.
#' When TRUE, the actual p-values are given.
#' @param seed Sets the seed in position_jitter(). Defaults to 1.
#' @param print.plot Boolean. Whether or not to print the plot(s). Defaults to FALSE.
#' @param add.title Boolean. Whether or not to add the name of the genomic feature to the plot. Defaults to FALSE.
#' @param wrap Boolean. Whether or not to wrap a list of plots into a single output. Defaults to FALSE.
#' @param add.legend Boolean. Whether or not to include a legend on the plots. Defaults to FALSE.
#' @param plot.titles A character vector containing the names of plot titles. Only used with stats.df is NULL and add.titles is TRUE.
#' 
#' @return A ggplot2 object.
#'
#' @import ggplot2
#' @import ggsignif
#' @import patchwork
#'
#' @export
#'
plot_global_overlap <- function(xint.obj, stats.df = NULL, numeric.p = FALSE, seed = 1, print.plot = FALSE,
                                add.title = FALSE, wrap = FALSE, add.legend = FALSE, plot.titles = NULL, 
                                y.scale = c("percent", "fraction"), ...) {

  if(!is.null(stats.df)) {
    if(!all(c("comparison", "p.adj") %in% colnames(stats.df))) {
      stop("stats.df must have columns 'comparison' and 'p.adj'.",
          call. = FALSE)
    }
  }

  if(is.list(xint.obj) && !is.data.frame(xint.obj)) {
    keep <- which(
      vapply(
        X = xint.obj,
        FUN = function(o){is(o, "xIntOverlap")},
        FUN.VALUE = logical(1)
        ))

    xint.obj <- xint.obj[keep]

    if(is.null(stats.df) && !is.null(plot.titles)) {
      if(length(plot.titles) != length(xint.obj)) {
        stop("The number of plot titles must be equal to the number of genomic features.",
             call. = FALSE)
      }
    }

    y.scale <- match.arg(y.scale)

    plots <- lapply(
      X = seq_along(xint.obj),
      FUN = function(i) {
        xo <- xint.obj[[i]]
        feat <- names(xint.obj)[i]
        if(!is.null(stats.df)) {
          tmp.stats <- stats.df[which(stats.df$feature == feat),]
        } else {
          tmp.stats <- NULL
        }
        base_overlap_plot(
          xint.obj = xo,
          stats.df = tmp.stats,
          numeric.p = numeric.p,
          seed = seed,
          print.plot = print.plot,
          add.title = add.title,
          add.legend = add.legend,
          plot.titles = plot.titles[i],
          y.scale = y.scale
        )
      }
    )

    names(plots) <- unique(stats.df$feature)

    if(wrap) {
      plots <- wrap_plots(plots, ...)
    }

  } else{
    if(is.null(stats.df) && !is.null(plot.titles)) {
      if(length(plot.titles) != 1) {
        stop("The number of plot titles must be equal to the number of genomic features.",
             call. = FALSE)
      }
    }

    plots <- base_overlap_plot(
      xint.obj = xint.obj,
      stats.df = stats.df,
      numeric.p = numeric.p,
      seed = seed,
      print.plot = print.plot,
      add.title = add.title,
      add.legend = add.legend,
      plot.titles = plot.titles,
      y.scale = y.scale
    )
  }

  if(print.plot) {
    print(plots)
  }

  return(plots)
}

base_overlap_plot <- function(xint.obj, stats.df = NULL, numeric.p = FALSE, seed = 1, print.plot = FALSE,
                              add.title = FALSE, add.legend = TRUE, plot.titles = NULL, y.scale = "percent") {
  plot.dat <- colData(xint.obj)

  if(y.scale == "percent") {
    plot.dat$fraction.overlap <- plot.dat$fraction.overlap * 100
  }

  if(!is.null(stats.df)) {
    stats.df$anno <- ifelse(stats.df$p.adj >= 0.05, "ns", 
                        ifelse(stats.df$p.adj < 0.05 & stats.df$p.adj >= 0.01, "*",
                                ifelse(stats.df$p.adj < 0.01 & stats.df$p.adj >= 0.001, "**",
                                      ifelse(stats.df$p.adj < 0.001, "***", NA))))
    
    annotation.pos <- calculate_annotation_positions(
      dat = plot.dat,
      stats.df = stats.df,
      y.scale = y.scale,
      min.gap = ifelse(y.scale == "percent", 2, 0.15)
    )

    plot.limits <- calculate_plot_limits(
      dat = plot.dat, 
      stats.df = stats.df, 
      annotation.positions = annotation.pos
    )
    y.max <- plot.limits[2]
  } else {
    y.max <- ifelse(y.scale == "percent", 100, 1)
  }

  p <- ggplot(plot.dat) +
    aes(x = condition, y = fraction.overlap, fill = condition) +
    geom_bar(stat = "summary", fun = "mean", color = "black", alpha = 0.8) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                 geom = "errorbar", width = 0.1, color = "black") +
    geom_point(
      shape = 21,
      size = 2.5,
      position = position_jitter(width = 0.2, seed = seed), 
      color = "black",
      show.legend = FALSE
      )
    
    if(!is.null(stats.df)) {
      p <- p + 
        ggsignif::geom_signif(
          comparisons = strsplit(stats.df$comparison, "-"),
          map_signif_level = FALSE,
          y_position = annotation.pos,
          annotations = if(numeric.p) {
                          formatC(stats.df$p.adj, format = "e", digits = 2)
                        } else {
                          stats.df$anno
                        },
          textsize = 4,
          step_increase = 0
        )
    }
    if(y.scale == "percent"){
      y.breaks <- seq(from = 0, to = 100, by = 20)
      base.max <- 100
    } else {
      y.breaks = seq(from = 0, to = 1, by = 0.2)
      base.max <- 1
    }
  
    y.label <- ifelse(y.scale == "percent", "Percent", "Fraction")
    
    p <- p +
      scale_fill_manual(values = get_custom_palette(length(unique(plot.dat$condition)))) +
      theme_bw() + 
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.key = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = ifelse(isTRUE(add.legend), "top", "none"),
        legend.box.margin = margin(-10,-10,-10,-10),
        legend.text = element_text(size = 10),
        legend.title = element_blank()
      ) +
      scale_y_continuous(limits = c(0, max(y.max + 0.05, base.max)), breaks = y.breaks) +
      labs(x = "Condition", y = paste0(y.label, " Overlap"))
  
  if(add.title) {
    if(!is.null(stats.df)) {
      feat <- unique(stats.df$feature)
        if(length(feat) > 1){
          stop("Statistics for multiple features are present in a single plot.", call. = FALSE)
        }
      } else {
      feat <- plot.titles
    }
    p <- p + labs(title = feat)
  }

  invisible(p)
}

calculate_annotation_positions <- function(dat, stats.df, y.scale = "percent", 
                                           base.offset = 0.05, step.size = 0.1, 
                                           min.gap = 2) {
  data.max <- max(dat$fraction.overlap)

  if(y.scale == "percent") {
    scale.range <- 100
  } else {
    scale.range <- 1
  }
  
  # Calculate base position (percentage of scale range above max data)
  base.position <- data.max + (scale.range * base.offset)
  
  # Calculate positions with adaptive spacing
  n.comparisons <- nrow(stats.df)
  step.increment <- scale.range * step.size
  
  positions <- base.position + seq(0, (n.comparisons - 1) * step.increment, 
                                   length.out = n.comparisons)
  
  # Ensure minimum gap between annotations
  if (length(positions) > 1) {
    for (i in 2:length(positions)) {
      min.pos <- positions[i-1] + min.gap
      if (positions[i] < min.pos) {
        positions[i] <- min.pos
      }
    }
  }
  
  return(positions)
}

calculate_plot_limits <- function(dat, stats.df, annotation.positions, 
                                  buffer = 0.05) {
  
  base.limit <- ifelse(max(dat$fraction.overlap) <= 1, 1, 100)
  
  if (!is.null(stats.df) && length(annotation.positions) > 0) {
    annotation.max <- max(annotation.positions)
    y.max <- annotation.max + (base.limit * buffer)
  } else {
    y.max <- base.limit
  }
  
  return(c(min = 0, max = y.max))
}
