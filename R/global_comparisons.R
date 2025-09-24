#' Wrap Basic Global Comparisons
#'
#' A basic wrapper for global comparison functions. Outputs a data frame containing effect size metrics and adjusted
#' p-values. Comparisons made via GLM utilize the quasi-Poisson model family and other default parameters in overlap_comparisons().
#' This function does not output all of the information returned in the comparison functions themselves. 
#' For more detailed output, run the functions individually.
#'
#' @param xint.obj An xIntOverlap object or a list of xIntOverlap objects (output by bulk_comparisons()).
#' @param comparison.type One of 'glm', 'mean', 'wilcoxon', or 'prop'. 
#' Determines the use of overlap_comparisons(), overlap_mean_tests(), wilcox_comparisons(), or overlap_prop_tests().
#' @param mean.test One of 't.test', or 'anova'. Only used when comparison.type = 'mean'.
#' @param prop.test One of 'fisher', 'chi', or 'G'. Only used when comparison.type = 'prop'.
#' @param p.adjust.method The method to use for p-value adjustments. Defaults to "BH".
#'
#' @return A data frame containing effect size metrics and p-values for each comparison.
#'
#' @examples
#' data(xobj)
#' # add examples
#'
#' @export
#'
global_comparisons <- function(xint.obj, comparison.type = c("glm", "mean", "wilcoxon", "prop"), 
                               mean.test = c("t.test", "anova"), prop.test = c("fisher", "chi", "G"), 
                               p.adjust.method = "BH") {
  comparison.type <- match.arg(comparison.type)

  if(comparison.type == "mean") {
    mean.test <- match.arg(mean.test)
  } else if(comparison.type == "prop") {
    prop.test = match.arg(prop.test)
  }

  p.adjust.method <- match.arg(p.adjust.method)

  if(is.list(xint.obj) && !is.data.frame(xint.obj)) {
    keep <- which(
      vapply(
        X = xint.obj,
        FUN = function(o){is(o, "xIntOverlap")},
        FUN.VALUE = logical(1)
        ))

    xint.obj <- xint.obj[keep]

    result.list <- lapply(
      X = xint.obj,
      FUN = function(xo) {
        if(comparison.type == "glm") {
          comparisons <- overlap_comparisons(xint.obj = xo, type = "global")
          comparisons <- contrast_stats(x = comparisons, contrast = Inf, p.adjust.method = p.adjust.method)
        } else if(comparison.type == "mean") {
          comparisons <- overlap_mean_tests(xint.obj = xo, method = mean.test, p.adjust.method = p.adjust.method)
        } else if(comparison.type == "wilcoxon") {
          comparison <- wilcox_comparisons(xint.obj = xo, type = "global")
          comparisons <- contrast_stats(x = comparisons, contrast = Inf, p.adjust.method = p.adjust.method)
        } else if(comparison.type == "prop") {
          comparisons <- overlap_prop_tests(xint.obj = xo, comparison = "pooled", p.adjust.method = p.adjust.method)
        }
        
        if(!"p.val" %in% colnames(comparisons)) {
          comparisons$p.val <- comparisons$p.adj
        }

        comparisons <- comparisons[ , c("comparison", "p.val", "p.adj")]
        roc <- overlap_effect(xint.obj = xo)
        tmp <- merge(roc, comparisons, by = "comparison", all = TRUE)
        return(tmp)
      }
    )

    result <- do.call(rbind, result.list)
  } else{
    if(comparison.type == "glm") {
      comparisons <- overlap_comparisons(xint.obj = xint.obj, type = "global")
      comparisons <- contrast_stats(x = comparisons, contrast = Inf, p.adjust.method = p.adjust.method)
    } else if(comparison.type == "mean") {
      comparisons <- overlap_mean_tests(xint.obj = xint.obj, method = mean.test, p.adjust.method = p.adjust.method)
    } else if(comparison.type == "wilcoxon") {
      comparisons <- wilcox_comparisons(xint.obj = xint.obj, type = "global")
      comparisons <- contrast_stats(x = comparisons, contrast = Inf, p.adjust.method = p.adjust.method)
    } else if(comparison.type == "prop") {
      comparisons <- overlap_prop_tests(xint.obj = xo, comparison = "pooled", p.adjust.method = p.adjust.method)
    }

    if(!"p.val" %in% colnames(comparisons)) {
      comparisons$p.val <- comparisons$p.adj
    }

    comparisons <- comparisons[ , c("comparison", "p.val", "p.adj")]

    roc <- overlap_effect(xint.obj = xint.obj)

    result <- merge(roc, comparisons, by = "comparison", all = TRUE)
  }

  result$feature <- gsub(".\\d+", "", row.names(result))
  return(result)
}