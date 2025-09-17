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
#' 
#' @return A ggplot2 object.
#'
#' @import ggplot2
#' @import ggsignif
#'
#' @export
#'
plot_global_overlap <- function(xint.obj, stats.df = NULL, numeric.p = FALSE, seed = 1) {
  if(!validObject(xint.obj)) {
    stop("xint.obj is not a valid xIntObject.", call. = FALSE)
  }

  if(!is.null(stats.df)) {
    if(!all(c("comparison", "p.adj") %in% colnames(stats.df))) {
      stop("stats.df must have columns 'comparison' and 'p.adj'.",
          call. = FALSE)
    }
  }

  plot.dat <- colData(xint.obj)

  if(!is.null(stats.df)) {
    stats.df$anno <- ifelse(stats.df$p.adj >= 0.05, "ns", 
                        ifelse(stats.df$p.adj < 0.05 & stats.df$p.adj >= 0.01, "*",
                                ifelse(stats.df$p.adj < 0.01 & stats.df$p.adj >= 0.001, "**",
                                      ifelse(stats.df$p.adj < 0.001, "***", NA))))
    
    y.max <- max(1, max(max(plot.dat$fraction.overlap) + 
                        seq(0.05, 0.06 * nrow(stats.df), length.out = nrow(stats.df))))
  } else {
    y.max <- 0.95
  }

  p <- ggplot(plot.dat) +
    aes(x = condition, y = fraction.overlap, fill = condition) +
    geom_bar(stat = "summary", fun = "mean", color = "black") +
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
          y_position = max(plot.dat$fraction.overlap) + 
                        seq(0.05, 0.06 * nrow(stats.df), length.out = nrow(stats.df)),
          annotations = if(numeric.p) {
                          formatC(stats.df$p.adj, format = "e", digits = 2)
                        } else {
                          stats.df$anno
                        },
          textsize = 4
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
      scale_y_continuous(limits = c(0, y.max + 0.05), breaks = seq(from = 0, to = 1, by = 0.2)) +
      labs(x = "Condition", y = "Fraction Overlap")

  print(p)  
  invisible(p)
}
  



