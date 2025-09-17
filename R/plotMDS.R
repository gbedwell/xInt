#' Make MDS Plot
#'
#' Create a multidimensional scaling (MDS) plot from a data matrix to visualize 
#' sample relationships based on leading fold changes or principal components.
#' The function calculates distances between samples and performs classical MDS
#' to reduce dimensionality for visualization.
#'
#' @param x A numeric matrix with features in rows and samples in columns, or an object that can be coerced to a matrix.
#' @param top Integer specifying the number of top variable features to use for distance calculations. Default is 500.
#' @param labels Character vector of sample labels. If NULL, column names of x are used, or sample indices if no column names exist.
#' @param dim.plot Numeric vector of length 2 specifying which dimensions to plot. Default is c(1,2).
#' @param n.dim Integer specifying the number of dimensions to calculate in MDS. Default is max(dim.plot).
#' @param feat.selection Character string specifying feature selection method. Either "pairwise" (top features selected for each sample pair) or "common" (same top features used for all comparisons). Default is "pairwise".
#' @param plot Logical indicating whether to print the plot. Default is TRUE.
#' @param dat Optional data frame containing sample metadata. If it contains a "condition" column, this will be used for coloring points.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot object containing the MDS plot. The plot is returned invisibly.
#'
plotMDS <- function(x, top = 500, labels = NULL, dim.plot = c(1,2), 
                    n.dim = max(dim.plot), feat.selection = "pairwise", 
                    plot = TRUE, dat = NULL, ...) {
  
  # Process input data
  x <- as.matrix(x)
  n.samples <- ncol(x)
  if(n.samples < 3) stop("Need at least 3 samples for MDS plot")
  
  # Remove rows with missing or Inf values
  bad <- rowSums(is.finite(x)) < n.samples
  if(any(bad)) x <- x[!bad, , drop = FALSE]
  n.probes <- nrow(x)
  
  # Use sample names as labels if not provided
  if(is.null(labels)) {
    labels <- colnames(x)
    if(is.null(labels)) labels <- 1:n.samples
  }
  
  # Validate parameters
  top <- min(top, n.probes)
  feat.selection <- match.arg(feat.selection, c("pairwise", "common"))
  
  # Calculate distance matrix
  dd <- matrix(0, nrow = n.samples, ncol = n.samples, dimnames = list(colnames(x), colnames(x)))
  
  if(feat.selection == "pairwise") {
    # Use top variable genes for each pair
    topindex <- n.probes - top + 1L
    for (i in 2:n.samples)
      for (j in 1:(i-1))
        dd[i,j] = sqrt(mean(sort.int((x[,i] - x[,j])^2, partial = topindex)[topindex:n.probes]))
    axislabel <- "Leading logFC dim"
  } else {
    # Use same top variable genes for all comparisons
    if(n.probes > top) {
      s <- rowMeans((x - rowMeans(x))^2)
      o <- order(s, decreasing=TRUE)
      x <- x[o[1:top], , drop=FALSE]
    }
    for (i in 2:n.samples)
      dd[i, 1:(i-1)] = sqrt(colMeans((x[,i] - x[,1:(i-1), drop=FALSE])^2))
    axislabel <- "Principal Component"
  }
  
  # Perform MDS
  mds_result <- suppressWarnings(cmdscale(as.dist(dd), k = n.dim))
  
  # Create plot data
  plot.dat <- data.frame(
    x = mds_result[, dim.plot[1]],
    y = mds_result[, dim.plot[2]],
    label = labels
  )
  
  # Add condition information from dat if provided
  if(!is.null(dat) && "condition" %in% colnames(dat)) {
    # Match conditions to samples
    if(nrow(dat) != n.samples) {
      warning("Number of rows in dat doesn't match number of samples. Using sample indices instead.")
      plot.dat$condition <- factor(1:n.samples)
    } else {
      plot.dat$condition <- factor(dat$condition)
    }
  } else {
    # If no condition provided, use sample indices
    plot.dat$condition <- factor(1:n.samples)
  }

  # Create ggplot
  p <- ggplot(plot.dat, aes(x = x, y = y, color = label, label = label)) +
      geom_text() +
      scale_color_manual(
        values = get_custom_palette(length(unique(plot.dat$condition)))
        ) +
      labs(
        x = paste(axislabel, dim.plot[1]),
        y = paste(axislabel, dim.plot[2]),
        title = "MDS Plot"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  if(plot) {
    print(p)
  }
  
  # Return plot object invisibly
  invisible(p)
}