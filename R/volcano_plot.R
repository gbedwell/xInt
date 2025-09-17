#' Volcano Plot for Per-feature Comparisons
#'
#' Creates a volcano plot showing -log10(p) vs log2 fold change or percent change versus for a specific contrast.
#' Takes as input the output of contrast_stats().
#'
#' @param x Output of contrast_stats()
#' @param p.cutoff p-value cutoff for highlighting significant features
#' @param fc.cutoff Log fold change cutoff for highlighting features with large effects
#' @param pc.cutoff Percent change cutoff for highlight features with large effects
#' @param p.col The column corresponding to the p-values of interest.
#' One of 'p.val' or 'p.adj' (default)
#' @param p.col The column corresponding to the effect size value of interest.
#' One of 'log2.fc' (default) or 'perc.change'
#' @param label.top Number of top features to label in the plot.
#' @param x.lim x-axis limits
#' @param y.lim y-axis limits
#'
#' @return A ggplot object or a list of ggplot objects (when multiple contrasts are present).
#'
#' @examples
#' data(xobj)
#' # add examples
#'
#' @import ggplot2
#' @import ggrepel
#' @export
#'
volcano_plot <- function(x, contrast, p.cutoff = 0.05, effect.cutoff = 1, p.col = "p.adj", 
                        effect.col = "log2.fc", label.top = 0, x.lim = c(-5,5), y.lim = c(0,6)) {

  if(!p.col %in% c("p.val", "p.adj")) {
    stop("p.col must be one of 'p.val' or 'p.adj'.", call. = FALSE)
  } else{
    if(p.col == "p.adj") {
      p.label = "FDR"
    } else {
      p.label = "p-value"
    }
  }

  if(!effect.col %in% c("log2.fc", "perc.change")) {
    stop("effect.col must be one of 'log2.fc' or 'perc.change'.", call. = FALSE)
  } else {
    if(effect.col == "log2.fc") {
      effect.label = "fold-change"
    } else {
      effect.label = "percent change"
    }
  }

  if(length(unique(x$comparison)) > 1) {
    contrasts <- unique(x$comparison)
    plot.list <- lapply(
      X = contrasts,
      FUN = function(contrast.name) {
        dat <- x[x$comparison == contrast.name, ]
        pp <- create_volcano_plot(
          x = dat,
          p.cutoff = p.cutoff,
          effect.cutoff = effect.cutoff,
          p.col = p.col,
          effect.col = effect.col,
          p.label = p.label,
          label.top = label.top,
          x.lim = x.lim,
          y.lim = y.lim
        )
        plot + ggtitle(contrast.name)
      }
    )
    names(plot.list) <- contrasts
    return(plot.list)
  } else {
    pp <- create_volcano_plot(
      x = x, 
      p.cutoff = p.cutoff, 
      effect.cutoff = effect.cutoff, 
      p.col = p.col, 
      effect.col = effect.col, 
      effect.label = effect.label,
      p.label = p.label,
      label.top = label.top, 
      x.lim = x.lim, 
      y.lim = y.lim
    )
    return(pp)
  }
}

create_volcano_plot <- function(x, p.cutoff, effect.cutoff, p.col, effect.col, 
                                effect.label, p.label, label.top, x.lim, y.lim) {
  
  p.dat <- data.frame(
    feature = x$feature,
    effect = x[[effect.col]],
    p = -log10(x[[p.col]]),
    type = ifelse(x[[p.col]] < p.cutoff & x[[effect.col]] > effect.cutoff, "Up",
           ifelse(x[[p.col]] < p.cutoff & x[[effect.col]] < -effect.cutoff, "Down", "None"))
  )

  p.dat$type <- factor(p.dat$type, levels = c("Up", "Down", "None"))

  p.dat$to.label <- FALSE
  if(label.top > 0) {
    sig.indices <- which(p.dat$type != "None")
    if(length(sig.indices) > 0) {
      top.indices <- sig.indices[order(x[[p.col]][sig.indices])[1:min(label.top, length(sig.indices))]]
      p.dat$to.label[top.indices] <- TRUE
    }
  }

  # Create the volcano plot
  p <- ggplot(p.dat) + 
        aes(x = effect, y = p, color = type) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = -log10(p.cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-effect.cutoff, effect.cutoff), linetype = "dashed") +
    scale_color_manual(values = c("Up" = "darkred", "Down" = "darkblue", "None" = "gray50")) +
    labs(
      x = paste0("log2(", effect.label, ")"),
      y = paste0("-log10(", p.label, ")")
    ) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.position = "top",
          legend.box.margin = margin(-10,-10,-10,-10),
          legend.text = element_text(size = 10),
          legend.title = element_blank()) +
    scale_x_continuous(limits = x.lim) + 
    scale_y_continuous(limits = y.lim) +
    geom_text_repel(
      data = subset(p.dat, to.label == TRUE),
      aes(label = feature), 
      max.overlaps = 20, 
      box.padding = 0.5,
      force = 5,
      show.legend = FALSE
    )
  
  return(p)
}