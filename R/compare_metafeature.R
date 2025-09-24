#' Compare Meta-Feature Distributions
#'
#' Compares meta-feature distributions between conditions. 
#' Comparisons are made using permutation tests on the difference between mean or median values.
#'
#' @param x The output of metafeature() with 'collapse = TRUE'.
#' @param n.perm The number of permutations. Defaults to 1000.
#' @param sub.samp The number of positions to sub-sample. Defaults to 10000.
#' Sub-sampling is ignored when NULL.
#' @param seed The random seed for reproducibility. Defaults to 1.
#' @param stat The test statistic to compare. One of 'quantiles' (default), 'median', or 'mean'.
#' @param quantiles The quantiles to calculate when stat = 'quantiles'. Defaults to c(0.25, 0.5, 0.75).
#' @param p.adj.method The method used to adjust p-values for multiple comparisons. Defaults to 'BH'.
#' @param bins The number of bins to use for binning the data. Defaults to 100. When NULL, operates on raw values.
#' @param combine.method Method for combining p-values when using multiple quantiles. One of 'stouffer' or 'fisher' (default).
#' 
#' @return A data frame holding the observed difference, p-value, and adjusted p-value for each comparison.
#'
#' @importFrom utils combn
#' @importFrom stats p.adjust median quantile pchisq pnorm qnorm
#'
#' @export
#'

compare_metafeature <- function(x, n.perm = 1000, sub.samp = 10000, seed = 1, 
                                stat = c("quantiles", "median", "mean"), 
                                quantiles = c(0.25, 0.5, 0.75), p.adj.method = "BH",
                                bins = NULL, combine.method = c("fisher", "stouffer")) {

  if(is.list(x) && !is.data.frame(x)) {
    stop("Input must be a data frame. Re-run metafeature() with 'collapse = TRUE'.", call. = FALSE)
  }
  if(!all(c("rel.position", "dataset", "condition") %in% colnames(x))) {
    stop("Input does not appear to be the output of metadata().", call. = FALSE)
  }
  
  stat <- match.arg(stat)
  combine.method <- match.arg(combine.method)
  
  if(!is.null(bins)) {
    if(!is.numeric(bins) || bins <= 0) {
      stop("bins must be a positive numeric value or NULL.", call. = FALSE)
    }
  }

  if(stat == "quantiles") {
    if(!is.numeric(quantiles) || any(quantiles < 0) || any(quantiles > 1)) {
      stop("Quantiles must be numeric values between 0 and 1.", call. = FALSE)
    }
  }
  
  # Function to combine p-values
  combine_pvalues <- function(p.vals, method = combine.method) {
    if(method == "fisher") {
      # Fisher's method
      # Handle p-values of 0 by setting them to a small value
      p.vals[p.vals == 0] <- .Machine$double.eps
      test.stat <- -2 * sum(log(p.vals))
      df <- 2 * length(p.vals)
      combined.p <- pchisq(test.stat, df, lower.tail = FALSE)
    } else if(method == "stouffer") {
      # Stouffer's method
      # Convert p-values to Z-scores
      z.scores <- qnorm(1 - p.vals)
      # Calculate combined Z-score
      combined.z <- sum(z.scores) / sqrt(length(z.scores))
      # Convert back to p-value
      combined.p <- 1 - pnorm(combined.z)
    }
    return(combined.p)
  }
  
  set.seed(seed)

  all.conditions <- unique(x$condition)
  condition.pairs <- combn(all.conditions, 2, simplify = FALSE)
  
  # Create bins and convert to bin centers for analysis
  bin_data <- function(values, n.bins) {
    breaks <- seq(0, 1, length.out = n.bins + 1)
    bin.centers <- (breaks[-1] + breaks[-length(breaks)]) / 2
    bins <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    # Return the bin centers for each value
    return(bin.centers[bins])
  }

  calc_stat <- function(values, method) {
    if(method == "median") {
      return(median(values))
    } else if(method == "mean") {
      return(mean(values))
    } else if(method == "quantiles") {
      return(quantile(values, probs = quantiles))
    }
  }

  calc_diff <- function(stat1, stat2) {
    if(is.vector(stat1) && length(stat1) > 1) {
      # For quantiles, return the differences for each quantile
      return(stat2 - stat1)
    } else {
      # For mean/median, return simple difference
      return(stat2 - stat1)
    }
  }

  results <- lapply(
    X = condition.pairs,
    FUN = function(pair) {
      cond1 <- pair[1]
      cond2 <- pair[2]
      
      cond1.dat <- x[x$condition == cond1,]
      cond1.pos <- cond1.dat$rel.position
      if(!is.null(sub.samp)) {
        cond1.pos <- sample(cond1.pos, size = min(sub.samp, length(cond1.pos)), replace = FALSE)
      }
      
      cond2.dat <- x[x$condition == cond2,]
      cond2.pos <- cond2.dat$rel.position
      if(!is.null(sub.samp)) {
        cond2.pos <- sample(cond2.pos, size = min(sub.samp, length(cond2.pos)), replace = FALSE)
      }

      # Apply binning if requested
      if(!is.null(bins)) {
        cond1.pos <- bin_data(cond1.pos, bins)
        cond2.pos <- bin_data(cond2.pos, bins)
      }
      
      # Calculate statistics
      cond1.stat <- calc_stat(cond1.pos, stat)
      cond2.stat <- calc_stat(cond2.pos, stat)
      
      # Calculate observed difference
      obs.diff <- calc_diff(cond1.stat, cond2.stat)
      
      # Combine data for permutation test
      all.pos <- c(cond1.pos, cond2.pos)
      n.cond1 <- length(cond1.pos)
      
      if(stat == "quantiles") {
        # For quantiles, calculate separate p-values for each quantile
        quantile.pvals <- numeric(length(quantiles))
        quantile.diffs <- obs.diff
        names(quantile.diffs) <- paste0("q", quantiles)
        
        for(i in seq_along(quantiles)) {
          # Run permutation test for this quantile
          perm.diffs <- numeric(n.perm)
          
          for(j in 1:n.perm) {
            shuf.pos <- sample(all.pos, size = length(all.pos), replace = FALSE)
            perm.cond1 <- shuf.pos[1:n.cond1]
            perm.cond2 <- shuf.pos[(n.cond1 + 1):length(all.pos)]
            
            perm.q1 <- quantile(perm.cond1, probs = quantiles[i])
            perm.q2 <- quantile(perm.cond2, probs = quantiles[i])
            
            perm.diffs[j] <- perm.q2 - perm.q1
          }
          
          # Calculate p-value for this quantile
          quantile.pvals[i] <- sum(abs(perm.diffs) >= abs(obs.diff[i])) / n.perm
        }
        
        # Combine p-values
        p.val <- combine_pvalues(quantile.pvals, combine.method)
        
        # Create result with individual quantile differences
        result <- data.frame(
          comparison = paste0(cond1, "-", cond2),
          stringsAsFactors = FALSE
        )
        
        for(i in seq_along(quantiles)) {
          result[[paste0("q", quantiles[i])]] <- quantile.diffs[i]
        }

        result$p.value <- p.val
        
      } else {
        # For mean/median, run a single permutation test
        perm.diffs <- numeric(n.perm)
        
        for(j in 1:n.perm) {
          shuf.pos <- sample(all.pos, size = length(all.pos), replace = FALSE)
          perm.cond1 <- shuf.pos[1:n.cond1]
          perm.cond2 <- shuf.pos[(n.cond1 + 1):length(all.pos)]
          
          if(stat == "mean") {
            perm.stat1 <- mean(perm.cond1)
            perm.stat2 <- mean(perm.cond2)
          } else {
            perm.stat1 <- median(perm.cond1)
            perm.stat2 <- median(perm.cond2)
          }
          
          perm.diffs[j] <- perm.stat2 - perm.stat1
        }
        
        p.val <- sum(abs(perm.diffs) >= abs(obs.diff)) / n.perm
        
        result <- data.frame(
          comparison = paste0(cond1, "-", cond2),
          difference = obs.diff,
          p.value = p.val,
          stringsAsFactors = FALSE
        )
      }

      return(result)
    }
  )

  results <- do.call(rbind, results)
  results$p.adj <- p.adjust(results$p.value, p.adj.method)
  rownames(results) <- results$comparison

  return(results)
}