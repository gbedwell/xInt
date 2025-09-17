#' Feature Overlap Comparisons for Single Samples
#'
#' Perform per-feature or feature-set integration targeting comparisons using Poisson exact tests.
#' This function is designed for cases without replicates per condition.
#' Library size offsets are used to account for differences in sequencing depth.
#'
#' @param xint.obj The xInt object of interest.
#' @param min.count Used to filter out lowly targeted features.
#'   If not provided, defaults to 1.
#' @param min.total.count Used alongside min.count to filter lowly targeted features.
#'   If not provided, defaults to 75% of the number of samples.
#' @param remove.zeros Boolean. Whether or not to remove rows with all zeros from the count matrix.
#'   Defaults to TRUE.
#' @param plot Boolean. Whether or not to output diagnostic plots.
#'   Defaults to TRUE.
#' @param contrasts Character vector of contrasts to test, in the form "A-B". If NULL (default),
#'   all pairwise comparisons between conditions will be created.
#' @param type Character string specifying the type of comparisons: per-feature ('local') or feature-set ('global').
#'   Default 'local'.
#' @param prior.count The constant added to count values to avoid taking log2(0).
#' @param lib.prior The constant added to library size when calculating log-CPM. 
#' When NULL, lib.prior = prior.count.
#'
#' @return A list containing the results of Poisson exact tests for each feature, including:
#'   \item{kept.feats}{Vector of features that passed filtering}
#'   \item{perc.change}{Matrix of percentage changes for pairwise comparisons}
#'   \item{log2.fc}{Matrix of log2 fold changes for pairwise comparisons}
#'   \item{p.values}{Matrix of p-values for pairwise comparisons}
#'   \item{contrasts}{Matrix of log fold changes for pairwise comparisons}
#'
#' @examples
#' data(xobj)
#' single_comp <- overlap_comparisons_single(xint.obj = xobj, plot = FALSE)
#'
#' @importFrom stats poisson.test
#' @import SummarizedExperiment
#'
#' @export
#'
single_overlap_comparisons <- function(xint.obj, min.count, min.total.count, 
                                      remove.zeros = TRUE, contrasts = NULL, 
                                      type = "local", prior.count = 0.5, lib.prior = NULL) {
  
  if(!validObject(xint.obj)) {
    stop("xint.obj is not a valid xIntObject.", call. = FALSE)
  }
  
  dat <- colData(xint.obj)

  if(length(unique(dat$condition)) != nrow(dat)){
    warning("Replication detected. For replicated data, consider overlap_comparisons().")
  }
  
  # Get count data based on type
  if(type == "local") {
    cs <- assay(xint.obj)
  } else {
    # For global analysis, use overlapping.sites from colData
    if(!"overlapping.sites" %in% colnames(dat)) {
      stop("For global analysis, 'overlapping.sites' column must exist in colData", call. = FALSE)
    }
    
    # Create a matrix with a single row for global analysis
    # This represents the total overlapping sites across all samples
    cs <- matrix(dat$overlapping.sites, nrow = 1, ncol = length(dat$overlapping.sites))
    rownames(cs) <- "global"
    colnames(cs) <- rownames(dat)
  }
  
  # Check if offset variable exists
  if(!"total.sites" %in% colnames(dat)) {
    stop("total.sites not found in colData.", call. = FALSE)
  }
  
  if(missing(min.count)) {
    min.count <- 1
  }
  
  if(missing(min.total.count)) {
    min.total.count <- floor(ncol(cs) * 0.75)
  }
  
  # Filter low count features
  keep <- rowSums(cs >= min.count) >= min.total.count
  filtered.cs <- cs[keep, , drop = FALSE]
  
  if(isTRUE(remove.zeros)) {
    # Remove rows with all zeros
    keep_nonzero <- rowSums(filtered.cs) > 0
    filtered.cs <- filtered.cs[keep_nonzero, , drop = FALSE]
  }
  
  # Get unique conditions
  conditions <- levels(dat$condition)
  if(length(conditions) < 2) {
    stop("At least two different conditions are required for comparisons.", call. = FALSE)
  }
  
  # Create contrast names
  if(is.null(contrasts)) {
    # Create all pairwise contrasts
    contrast.names <- character()
    for(i in seq(1, (length(conditions) - 1))) {
      for(j in seq((i + 1), length(conditions))) {
        if(conditions[i] != conditions[j]) {
          contrast.name <- paste(conditions[i], "-", conditions[j], sep = "")
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
  
  if(is.null(lib.prior)) {
    lib.prior = prior.count
  }

  lcpm <- calculate_cpm(
    filtered.cs, 
    lib.size = dat$total.sites, 
    log = TRUE, 
    prior.count = prior.count,
    lib.prior = lib.prior
    )

  contrast.mat.results <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  percent.change.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  log2.fc.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  p.values.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  
  rownames(contrast.mat.results) <- rownames(p.values.mat) <- 
    rownames(percent.change.mat) <- rownames(log2.fc.mat) <- rownames(filtered.cs)
  
  colnames(contrast.mat.results) <- colnames(p.values.mat) <- 
    colnames(percent.change.mat) <- colnames(log2.fc.mat) <- contrast.names
  
  # Process each contrast
  for(c in 1:n.contrasts) {
    # Parse contrast
    contrast_parts <- strsplit(contrast.names[c], "-")[[1]]
    cond1 <- contrast_parts[1]
    cond2 <- contrast_parts[2]
    
    # Get samples for each condition
    samples1 <- rownames(dat)[dat$condition == cond1]
    samples2 <- rownames(dat)[dat$condition == cond2]
    
    if(length(samples1) == 0 || length(samples2) == 0) {
      stop("No samples found for condition ", ifelse(length(samples1) == 0, cond1, cond2))
    }
    
    # Get library sizes (offsets)
    lib.size1 <- sum(dat[samples1, "total.sites"])
    lib.size2 <- sum(dat[samples2, "total.sites"])
    
    # For each feature, perform Poisson exact test
    for(i in 1:n.features) {
      # Get counts
      counts1 <- sum(filtered.cs[i, samples1, drop = FALSE])
      counts2 <- sum(filtered.cs[i, samples2, drop = FALSE])
      
      # Calculate normalized rates (per unit of library size)
      rate1 <- counts1 / lib.size1
      rate2 <- counts2 / lib.size2
      
      # Calculate log fold change with pseudocount for consistency
      rate1.adj <- (counts1 + prior.count) / (lib.size1 + lib.prior)
      rate2.adj <- (counts2 + prior.count) / (lib.size2 + lib.prior)
      log.fc <- log(rate2.adj / rate1.adj)
      
      # Store fold change results
      contrast.mat.results[i, c] <- log.fc
      percent.change.mat[i, c] <- (exp(log.fc) - 1) * 100
      log2.fc.mat[i, c] <- log.fc / log(2)
      
      # Perform Poisson exact test
      test_result <- stats::poisson.test(
        x = c(counts1, counts2),
        T = c(lib.size1, lib.size2),
        alternative = "two.sided",
        conf.level = 0.95
      )
      p.value <- test.result$p.value
      
      p.values.mat[i, c] <- p.value
    }
  }
  
  # Create a dummy contrast matrix to match the structure expected by contrast_hits
  contrast.mat <- matrix(0, nrow = length(contrast.names), ncol = length(contrast.names))
  rownames(contrast.mat) <- colnames(contrast.mat) <- contrast.names
  diag(contrast.mat) <- 1
  
  # Return results with the structure expected by contrast_hits
  results <- list(
    kept.feats = rownames(filtered.cs),
    log.cpm = lcpm,
    perc.change = percent.change.mat,
    log2.fc = log2.fc.mat,
    p.values = matrix(NA, nrow = n.features, ncol = 0),  # Empty matrix for coefficient p-values
    contrasts = contrast.mat.results,
    contrast.p.vals = p.values.mat,  # Rename p.values to contrast.p.vals
    model.info = list(
      formula = as.formula("~ condition"),
      family = "poisson.test",
      design = NULL
    ),
    contrast.mat = contrast.mat
  )
  
  return(results)
}