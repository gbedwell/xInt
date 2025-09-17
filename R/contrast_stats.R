#' Extract Comparison Statistics
#'
#' Neatly extract p-value and effect size information for specified contrasts from the
#' outputs of overlap_comparisons(), single_overlap_comparisons(), or wilcox_overlaps().
#'
#' @param x The output from overlap_comparisons() or single_overlap_comparisons()
#' @param contrast The name or numeric index of the contrast to extract.
#'   If Inf, results for all contrasts will be returned in a single dataframe.
#' @param n.hits Number of top hits to extract. Default is Inf (all hits).
#' @param p.cutoff P-value cutoff for filtering results. Default is 0.05.
#' @param fc.cutoff Log2 fold change cutoff for filtering results. Default is 1.
#'   Only applied when contrast is specified.
#' @param pc.cutoff Percent change cutoff for filtering results. Default is 0.
#'   Only applied when contrast is specified.
#' @param p.adjust.method Method for p-value adjustment. Default is "BH" (Benjamini-Hochberg).
#' @param sort.by Column to sort results by. Default is "p.adj".
#' @param decreasing Boolean indicating if sorting should be in decreasing order. Default is FALSE.
#' @param significant.only Boolean indicating whether or not to only return significant features. Defaults to FALSE.
#'
#' @return A data frame with rows for each feature and columns for various statistics
#'
#' @examples
#' data(xobj)
#' glm.results <- overlap_comparisons(xint.obj = xobj)
#' 
#' contrast_stats(
#'   glm.results,
#'   contrast = 1,
#'   n.hits = 2
#'   p.cutoff = 1
#'  )
#'
#' @export
#'
contrast_stats <- function(x, contrast = Inf, n.hits = Inf, 
                          p.cutoff = 0.05, fc.cutoff = 1, pc.cutoff = 0,
                          p.adjust.method = "BH", sort.by = "p.adj",
                          decreasing = FALSE, significant.only = FALSE) {

  if(!all(c("kept.feats", "log.cpm", "perc.change", "log2.fc", "p.values", "contrasts") %in% names(x))) {
    stop("Ensure the provided input is the output of overlap_comparisons() or single_overlap_comparisons().")
  }
  
  if(is.null(p.adjust.method)) {
    use.raw.p <- TRUE
  } else {
    use.raw.p <- FALSE
  }

  if(fc.cutoff > 0 && pc.cutoff > 0) {
    stop("Only one of fc.cutoff or pc.cutoff can be specified.", call. = FALSE)
  }
  
  # Working with contrasts
  contrast.names <- colnames(x$contrasts)
  
  if(length(contrast.names) == 0) {
    stop("No contrasts found in results")
  }

  if(is.infinite(contrast) || length(contrast) > 1) {
    if(!is.infinite(contrast)) {
      contrast.names <- contrast.names[contrast]
    }

    all.results <- lapply(
      X = seq_along(contrast.names),
      FUN = function(i) {
        results.df <- data.frame(
          feature = x$kept.feats,
          comparison = contrast.names[i],
          log2.fc = x$log2.fc[, i],
          perc.change = x$perc.change[, i],
          p.val = x$contrast.p.vals[, i],
          stringsAsFactors = FALSE
        )

        if(!is.null(p.adjust.method)) {
          results.df$p.adj <- p.adjust(results.df$p.val, method = p.adjust.method)
        } else {
          results.df$p.adj <- results.df$p.val
        }
        
        p.col <- ifelse(use.raw.p, "p.val", "p.adj")
        if(fc.cutoff == 0 && pc.cutoff == 0) {
          results.df$significant <- results.df[[p.col]] <= p.cutoff
        } else if(fc.cutoff > 0) {
          results.df$significant <- results.df[[p.col]] <= p.cutoff & abs(results.df$log2.fc) >= fc.cutoff
        } else if(pc.cutoff > 0) {
          results.df$significant <- results.df[[p.col]] <= p.cutoff & abs(results.df$perc.change) >= pc.cutoff
        }

        return(results.df)
      }
    )

    combined.results <- do.call(rbind, all.results)

    if(sort.by %in% c("p.adj", "p.val")) {
      combined.results <- combined.results[order(combined.results[[sort.by]], decreasing = decreasing), ]
    } else if(sort.by %in% c("log2.fc", "perc.change")) {
      combined.results <- combined.results[order(abs(combined.results[[sort.by]]), decreasing = decreasing), ]
    }

    if(isTRUE(significant.only)){
      combined.results <- combined.results[combined.results$significant == TRUE, ]
    }

    if(is.finite(n.hits) && n.hits < nrow(combined.results)) {
      combined.results <- combined.results[1:n.hits, ]
    }

    rownames(combined.results) <- NULL

    return(combined.results)
  }

  # Handle numeric index
  if(is.numeric(contrast)) {
    if(contrast < 1 || contrast > length(contrast.names)) {
      stop("Contrast index out of range. Must be between 1 and ", length(contrast.names))
    }
    idx <- contrast
    contrast.name <- contrast.names[idx]
  } else {
    # Check if contrast exists by name
    if(!contrast %in% contrast.names) {
      stop("Contrast '", contrast, "' not found in results. Available contrasts: ", 
            paste(contrast.names, collapse = ", "))
    }
    # Get contrast index
    idx <- which(contrast.names == contrast)
    contrast.name <- contrast
  }
  
  # Create data frame with results
  results.df <- data.frame(
    feature = x$kept.feats,
    comparison = contrast.name,
    log2.fc = x$log2.fc[, idx],
    perc.change = x$perc.change[, idx],
    p.val = x$contrast.p.vals[, idx],
    stringsAsFactors = FALSE
  )
  
  results.df$p.adj <- p.adjust(results.df$p.val, method = p.adjust.method)

  p.col <- ifelse(use.raw.p, "p.val", "p.adj")
  if(fc.cutoff == 0 && pc.cutoff == 0) {
    results.df$significant <- results.df[[p.col]] <= p.cutoff
  } else if(fc.cutoff > 0) {
    results.df$significant <- results.df[[p.col]] <= p.cutoff & abs(results.df$log2.fc) >= fc.cutoff
  } else if(pc.cutoff > 0) {
    results.df$significant <- results.df[[p.col]] <= p.cutoff & abs(results.df$perc.change) >= pc.cutoff
  }

  # Sort results
  if(sort.by %in% c("p.adj", "p.val")){
    results.df <- results.df[order(results.df[[sort.by]], decreasing = decreasing), ]
  } else if(sort.by %in% c("log2.fc", "perc.change")) {
    results.df <- results.df[order(abs(results.df[[sort.by]]), decreasing = decreasing), ]
  }

  if(isTRUE(significant.only)){
    results.df <- results.df[results.df$significant == TRUE, ]
  }

  # Limit to n.hits
  if(is.finite(n.hits) && n.hits < nrow(results.df)) {
    results.df <- results.df[1:n.hits, ]
  }

  rownames(results.df) <- NULL

  return(results.df)
}