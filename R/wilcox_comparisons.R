#' Feature Overlap Comparisons Using Wilcoxon
#'
#' Perform per-feature integration targeting comparisons using Wilcoxon rank-sum tests.
#' This is a non-parametric alternative to GLM-based approaches.
#' Uses CPM normalization to account for library size differences.
#'
#' @param xint.obj The xInt object of interest.
#' @param min.count Used to filter out lowly targeted features.
#'   If not provided, defaults to 1.
#' @param min.total.count Used alongside min.count to filter lowly targeted features.
#'   If not provided, defaults to one-half of the number of samples.
#' @param remove.zeros Boolean. Whether or not to remove rows with all zeros from the count matrix.
#'   Defaults to TRUE.
#' @param plot Boolean. Whether or not to output diagnostic plots.
#'   Defaults to TRUE.
#' @param contrasts Character vector of contrasts to test, in the form "A-B". If NULL (default),
#'   all pairwise comparisons between conditions will be created.
#' @param prior.count Prior count added for log-CPM calculation to avoid log(0).
#'   Defaults to 0.5 (more conservative than the GLM default of 2).
#'
#' @return A list containing the results of Wilcoxon tests for each feature, including:
#'   \item{kept.feats}{Vector of features that passed filtering}
#'   \item{log.cpm}{Log-CPM normalized expression matrix}
#'   \item{perc.change}{Matrix of percentage changes for pairwise comparisons}
#'   \item{log2.fc}{Matrix of log2 fold changes for pairwise comparisons}
#'   \item{contrasts}{Matrix of log fold changes for pairwise comparisons}
#'   \item{contrast.p.vals}{Matrix of p-values for pairwise comparisons}
#'   \item{model.info}{Information about the method used}
#'   \item{contrast.mat}{The contrast matrix used for comparisons}
#'
#' @examples
#' data(xobj)
#' wilcox_comparisons(xint.obj = xobj, plot = FALSE)
#'
#' @importFrom stats wilcox.test
#' @import SummarizedExperiment
#' @import ggplot2
#'
#' @export
#'
wilcox_comparisons <- function(xint.obj,
                               type = c("local", "global"),
                               min.count = 1, 
                               min.total.count = NULL,
                               remove.zeros = TRUE,
                               plot = TRUE,
                               contrasts = NULL) {
  
  if(!validObject(xint.obj)) {
    stop("xint.obj is not a valid xIntObject.", call. = FALSE)
  }

  type = match.arg(type)
  
  dat <- colData(xint.obj)
  
  if(type == "local") {
    cs <- assay(xint.obj)
  } else {
    if(!"overlapping.sites" %in% colnames(dat)) {
      stop("For global analysis, 'overlapping.sites' column must exist in colData", call. = FALSE)
    }
    
    cs <- matrix(dat$overlapping.sites, nrow = 1, ncol = length(dat$overlapping.sites))
    rownames(cs) <- "global"
    colnames(cs) <- rownames(dat)
  }
  
  if(!"total.sites" %in% colnames(dat)) {
    stop("total.sites not found in colData.", call. = FALSE)
  }

  if(missing(min.count)) {
    min.count <- 1
  }
  
  if(is.null(min.total.count)) {
    min.total.count <- floor(ncol(cs) * 0.75)
  }
  
  keep <- rowSums(cs >= min.count) >= min.total.count
  filtered.cs <- cs[keep, , drop = FALSE]
  
  if(isTRUE(remove.zeros)) {
    keep_nonzero <- rowSums(filtered.cs) > 0
    filtered.cs <- filtered.cs[keep_nonzero, , drop = FALSE]
  }
  
  prior.count <- 0.5
  lib.prior <- 1

  lcpm <- calculate_cpm(filtered.cs, lib.size = dat$total.sites, log = TRUE, 
                        prior.count = prior.count, lib.prior = lib.prior)
  
  # Generate diagnostic plots
  if(isTRUE(plot) && type == "local") {
    plotMDS(lcpm, labels = dat$condition, dat = dat)
    
    # For Wilcoxon, show distribution of log-CPM values
    plot_lcpm_distributions(lcpm, dat, prior.count)
  }
  
  conditions <- levels(dat$condition)
  if(length(conditions) < 2) {
    stop("At least two different conditions are required for comparisons.", call. = FALSE)
  }

  cond.lens <- vapply(
    X = conditions,
    FUN = function(x){
      length(dat$condition[dat$condition == x])
      },
    FUN.VALUE = numeric(1)
   )

  low.rep.conditions <- names(cond.lens[cond.lens < 3])
  if(length(low.rep.conditions) > 0) {
    stop("The following conditions have < 3 replicates: ",
          paste(low.rep.conditions, collapse = ", "), ".",
          call. = FALSE)
  }
  
  # Create contrast names
  if(is.null(contrasts)) {
    # Create all pairwise contrasts
    contrast.names <- character()
    for(i in seq(1, (length(conditions) - 1))) {
      for(j in seq((i + 1), length(conditions))) {
        if(conditions[i] != conditions[j]) {
          contrast.name <- paste(conditions[j], "-", conditions[i], sep = "")
          contrast.names <- c(contrast.names, contrast.name)
        }
      }
    }
  } else {
    contrast.names <- contrasts
  }
  
  # Initialize result matrices
  n.features <- nrow(filtered.cs)
  n.contrasts <- length(contrast.names)
  
  contrast.mat.results <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  percent.change.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  log2.fc.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  contrast.pval.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  
  rownames(contrast.mat.results) <- rownames(contrast.pval.mat) <- 
    rownames(percent.change.mat) <- rownames(log2.fc.mat) <- rownames(filtered.cs)
  
  colnames(contrast.mat.results) <- colnames(contrast.pval.mat) <- 
    colnames(percent.change.mat) <- colnames(log2.fc.mat) <- contrast.names
  
  # Process each contrast
  for(c in 1:n.contrasts) {
    # Parse contrast
    contrast.parts <- strsplit(contrast.names[c], "-")[[1]]
    cond1 <- contrast.parts[1]
    cond2 <- contrast.parts[2]
    
    # Get samples for each condition
    samples1 <- rownames(dat)[dat$condition == cond1]
    samples2 <- rownames(dat)[dat$condition == cond2]
    
    if(length(samples1) == 0 || length(samples2) == 0) {
      stop("No samples found for condition ", ifelse(length(samples1) == 0, cond1, cond2))
    }
    
    # For each feature, perform Wilcoxon test on log-CPM values
    for(i in 1:n.features) {
      # Get log-CPM values for each group
      group1.values <- lcpm[i, samples1]
      group2.values <- lcpm[i, samples2]
      
      # Calculate fold change from group medians
      # Using medians retains consistency with Wilcoxon
      median1 <- median(group1.values)
      median2 <- median(group2.values)
      log.fc <- median1 - median2  # log2(FC) = log2(A) - log2(B)
      
      # Store fold change results
      contrast.mat.results[i, c] <- NA
      log2.fc.mat[i, c] <- log.fc  # Already in log2 scale
      percent.change.mat[i, c] <- (2^log.fc - 1) * 100  # Convert from log2 to percentage
      
      # Perform Wilcoxon rank-sum test
      if(length(unique(c(group1.values, group2.values))) > 1) {
        test.result <- wilcox.test(
          x = group1.values,
          y = group2.values,
          alternative = "two.sided",
          exact = NULL
        )
        p.value <- test.result$p.value
      } else {
        p.value <- 1.0
      }
      contrast.pval.mat[i, c] <- p.value
    }
  }
  
  # Create contrast matrix for compatibility with other functions
  contrast.mat <- matrix(0, nrow = length(contrast.names), ncol = length(contrast.names))
  rownames(contrast.mat) <- colnames(contrast.mat) <- contrast.names
  diag(contrast.mat) <- 1
  
  # Return results with structure compatible with overlap_comparisons
  results <- list(
    kept.feats = rownames(filtered.cs),
    log.cpm = lcpm,
    perc.change = percent.change.mat,
    log2.fc = log2.fc.mat,
    coefficients = matrix(NA, nrow = n.features, ncol = 0),
    p.values = matrix(NA, nrow = n.features, ncol = 0),
    contrasts = contrast.mat.results,
    contrast.p.vals = contrast.pval.mat,
    model.info = list(
      method = "wilcoxon",
      normalization = "log2-CPM",
      prior.count = prior.count
    ),
    contrast.mat = contrast.mat
  )

  if(isTRUE(plot) && type == "local") {
    plot_pvals(ps = c(results$contrast.p.vals))
  }
  
  return(results)
}

#' Plot log-CPM distributions by condition
#'
#' @param lcpm Log-CPM matrix
#' @param dat Sample metadata
#' @param prior.count The prior.count added to values before log-transformation
#'
#'
plot_lcpm_distributions <- function(lcpm, dat, prior.count) {
  # Create long-format data for plotting
  plot.dat <- data.frame(
    sample = rep(colnames(lcpm), each = nrow(lcpm)),
    condition = rep(dat$condition, each = nrow(lcpm)),
    log.cpm = as.vector(lcpm)
  )
  
  # Remove infinite values for plotting
  plot.dat <- plot.dat[is.finite(plot.dat$log.cpm), ]
  
  p <- ggplot(plot.dat, aes(x = log.cpm, fill = condition)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(
        values = get_custom_palette(length(unique(plot.dat$condition)))
        ) +
    facet_wrap(~ condition, scales = "free_y") +
    labs(
      title = "log2-CPM Distributions",
      x = "log2(CPM)",
      y = "Density"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  print(p)
  invisible(p)
}

plot_pvals <- function(ps) {
  p <- ggplot(data = data.frame(p = ps), aes(x = p)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), color = "black", fill = "gray75") +
    labs(
      title = "Combined p-value Histogram",
      x = "p-value",
      y = "Count"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      axis.text = element_text(size=12),
      axis.title = element_text(size=14)
    )
  
  # Print the plot
  print(p)
  
  # Return the plot object invisibly
  invisible(p)
}