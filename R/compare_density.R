#' Compare Feature Density Between Conditions
#'
#' Performs statistical tests to compare feature density between different conditions.
#'
#' @param x The output of feature_density() with options average = TRUE and conditions are given.
#' @param p.adjust.method Method for p-value adjustment. Default is "BH" (Benjamini-Hochberg).
#' @param test Statistical test to use. Either 't.test' (default) or 'anova'.
#'
#' @return A data frame containing the results of statistical comparisons between conditions.
#'
#' @export
#'
compare_density <- function(x, p.adjust.method = "BH", test = c("t.test", "anova")) {
  
  test <- match.arg(test)
  
  if(!all(c("sample", "condition", "avg.density") %in% colnames(x))) {
    stop("x must match the output of feature_density().",
         call. = FALSE)
  }
  
  # Convert condition to factor if it's not already
  x$condition <- as.factor(x$condition)
  
  # Get unique conditions
  conditions <- levels(x$condition)

  # Count replicates per condition
  cond.lens <- vapply(
    X = conditions,
    FUN = function(c) {
      sum(x$condition == c)
    },
    FUN.VALUE = numeric(1)
  )
  
  multi.rep <- names(cond.lens[cond.lens > 1])
  single.rep <- names(cond.lens[cond.lens == 1])
  
  if(length(multi.rep) == 0) {
    stop("No replicate datasets are provided.",
         call. = FALSE)
  }
  
  if(test == "anova") {
    if(length(single.rep) == 0) {
      # All conditions have replicates, perform ANOVA
      dens.aov <- aov(avg.density ~ condition, data = x)
      aov.hsd <- TukeyHSD(dens.aov, conf.level = 0.95)
      comparisons <- data.frame(
        comparison = rownames(aov.hsd$condition),
        diff = aov.hsd$condition[, 1],
        lower = aov.hsd$condition[, 2],
        upper = aov.hsd$condition[, 3],
        p.val = aov.hsd$condition[, 4],
        p.adj = aov.hsd$condition[, 4]
      )
      
      rownames(comparisons) <- rownames(aov.hsd$condition)
    } else {
      stop("Cannot perform ANOVA if all data are not replicated.", call. = FALSE)
    }
  } else {
    # t-test comparisons
    comparisons <- list()
    
    if(length(multi.rep) > 1) {
      for(i in seq(1, (length(multi.rep) - 1))) {
        for(j in seq((i + 1), length(multi.rep))) {
          if(multi.rep[i] != multi.rep[j]) {
            s1 <- x$avg.density[x$condition == multi.rep[i]]
            s2 <- x$avg.density[x$condition == multi.rep[j]]
            contrast.name <- paste(multi.rep[i], "-", multi.rep[j], sep = "")
            tt <- t.test(x = s1,
                         y = s2,
                         alternative = "two.sided",
                         var.equal = TRUE)
            
            out <- data.frame(
              comparison = contrast.name,
              t = tt$statistic,
              df = tt$parameter,
              mu1 = tt$estimate[1],
              mu2 = tt$estimate[2],
              p.val = tt$p.value
            )
            
            comparisons[[contrast.name]] <- out
          }
        }
      }
    }
    
    if(length(single.rep) != 0) {
      for(i in seq(1, length(multi.rep))) {
        for(j in seq(1, length(single.rep))) {
          if(multi.rep[i] != single.rep[j]) {
            s1 <- x$avg.density[x$condition == multi.rep[i]]
            s2 <- x$avg.density[x$condition == single.rep[j]]
            contrast.name <- paste(multi.rep[i], "-", single.rep[j], sep = "")
            
            m1 <- mean(s1)
            n1 <- length(s1)
            m2 <- mean(s2)
            n2 <- length(s2)
            
            se <- sqrt((1/n1 + 1/n2) * ((n1 - 1) * sd(s1)^2 + (n2 - 1) * 0^2) / (n1 + n2 - 2))
            df <- n1 + n2 - 2
            t <- (m1 - m2) / se
            
            out <- data.frame(
              comparison = contrast.name,
              t = t,
              df = df,
              mu1 = m1,
              mu2 = m2,
              p.val = 2 * pt(-abs(t), df)
            )
            
            comparisons[[contrast.name]] <- out
          }
        }
      }
      
      if(length(single.rep) > 1) {
        warning("The following datasets cannot be compared: ",
                paste(single.rep, collapse = ", "), ".",
                call. = FALSE)
      }
    }
    
    comparisons <- do.call(rbind, comparisons)
    comparisons$p.adj <- p.adjust(p = comparisons$p.val, method = p.adjust.method)
  }
  
  return(comparisons)
}