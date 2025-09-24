#' Plot Density
#'
#' Plot per-condition average feature density.
#'
#' @param x The output of feature_density() with options average = TRUE and conditions are given.
#' @param stats.df The output of compare_density().
#' @param condition.levels The condition levels for the samples. If NULL (default), the levels are created in the order they appear.
#' @param print.plot Boolean. Whether or not to print the plot(s). Defaults to FALSE.
#' @param y.label The label of the y-axis. Defaults to "Genes / Mb".
#' @param numeric.p Boolean. When FALSE (default), significance is denoted symbolically.
#' p > 0.05 is 'ns', 0.01 >= p < 0.05 is '*', 0.001 >= p < 0.01 is '**', and p < 0.001 is '***'.
#' When TRUE, the actual p-values are given.
#' @param seed Sets the seed in position_jitter(). Defaults to 1.
#' 
#' @return A ggplot2 object.
#'
#' @import ggplot2
#'
#' @export
#'
plot_density <- function(x, stats.df = NULL, condition.levels = NULL, print.plot = FALSE, y.label = "Genes / Mb",
                         numeric.p = FALSE, seed = 1) {

  if(!all(c("sample", "condition", "avg.density") %in% colnames(x))) {
    stop("x must match the output of feature_density().",
         call. = FALSE)
  }

  if(!is.null(stats.df)) {
    if(!all(c("comparison", "p.adj") %in% colnames(stats.df))) {
      stop("stats.df must have columns 'comparison' and 'p.adj'.",
          call. = FALSE)
    }
  }

  if(!is.null(condition.levels)) {
    x$condition <- factor(x$condition, levels = condition.levels)
  } else {
    x$condition <- factor(x$condition, levels = unique(x$condition))
  }

  pp <- base_dens_plot(
    x = x,
    stats.df = stats.df, 
    print.plot = print.plot, 
    y.label = y.label,
    numeric.p = numeric.p,
    seed = seed
  )

}

base_dens_plot <- function(x, stats.df = NULL, numeric.p = FALSE, seed = 1, print.plot = FALSE, y.label = "Genes / Mb") {
  plot.dat <- x

  if(!is.null(stats.df)) {
    stats.df$anno <- ifelse(stats.df$p.adj >= 0.05, "ns", 
                        ifelse(stats.df$p.adj < 0.05 & stats.df$p.adj >= 0.01, "*",
                                ifelse(stats.df$p.adj < 0.01 & stats.df$p.adj >= 0.001, "**",
                                      ifelse(stats.df$p.adj < 0.001, "***", NA))))
    
    annotation.pos <- calculate_annotation_positions(
      dat = plot.dat,
      stats.df = stats.df
    )

    plot.limits <- calculate_plot_limits(
      dat = plot.dat, 
      stats.df = stats.df, 
      annotation.positions = annotation.pos
    )

    y.max <- plot.limits[2]

  } else {
    dens.max <- max(plot.dat$avg.density)
    y.max <- 10 * (dens.max %/% 10 + as.logical(dens.max %% 10))
  }

  p <- ggplot(plot.dat) +
    aes(x = condition, y = avg.density, fill = condition) +
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

    p <- p +
      scale_fill_manual(values = get_custom_palette(length(unique(plot.dat$condition)))) +
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
      scale_y_continuous(limits = c(0, y.max)) +
      labs(x = "Condition", y = y.label)

  if(print.plot) {
    print(p)
  }

  invisible(p)
}


calculate_annotation_positions <- function(dat, stats.df, base.offset = 0.05, 
                                           step.size = 0.1, min.gap = 2) {
  dens.max <- max(dat$avg.density)
  scale.range <- 10 * (dens.max %/% 10 + as.logical(dens.max %% 10))
  
  # Calculate base position (percentage of scale range above max data)
  base.position <- dens.max + (scale.range * base.offset)
  
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
  
  dens.max <- max(dat$avg.density)
  base.limit <- 10 * (dens.max %/% 10 + as.logical(dens.max %% 10))
  
  if (!is.null(stats.df) && length(annotation.positions) > 0) {
    annotation.max <- max(annotation.positions)
    y.max <- annotation.max + (base.limit * buffer)
  } else {
    y.max <- base.limit
  }
  
  return(c(min = 0, max = y.max))
}