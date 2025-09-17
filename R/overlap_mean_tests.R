#' Global Overlap Comparisons Using Mean Tests
#'
#' Perform global pairwise comparisons on an xIntOverlap object using tests for equal means.
#' Comparisons are made using either t-tests or ANOVA followed by Tukey HSD.
#'
#' @param xint.obj The xInt object of interest.
#' @param p.adjust.method The p-value adjustment method.
#' Defaults to "BH". See p.adjust() for more information.
#' @param method Method for comparison. One of "t.test" (default) for pairwise t-tests or "anova" for ANOVA followed by Tukey HSD.
#'
#' @return A data frame containing the results of each comparison.
#'
#' @examples
#' data(xobj)
#' mean_tests(xint.obj = xobj)
#'
#' @importFrom stats t.test p.adjust aov TukeyHSD wilcox.test
#' @import methods
#'
#' @export
#'
overlap_mean_tests <- function(xint.obj, p.adjust.method = "BH", method = c("t.test", "anova")){

  method <- match.arg(method)

  if(!validObject(xint.obj)){
    stop("xint.obj is not a valid xIntObject.",
         call. = FALSE)
    }

  dat <- colData(xint.obj)

  lcpm <- calculate_cpm(
    y = dat$overlapping.sites,
    lib.size = dat$total.sites,
    log = TRUE,
    prior.count = 0.5,
    lib.prior = 1
    )

  names(lcpm) <- dat$condition

  conditions <- levels(dat$condition)

  cond.lens <- vapply(
    X = conditions,
    FUN = function(x){
      length(dat$condition[dat$condition == x])
      },
    FUN.VALUE = numeric(1)
   )

  multi.rep <- names(cond.lens[cond.lens > 1])
  single.rep <- names(cond.lens[cond.lens == 1])

  if(length(multi.rep) == 0){
    stop("No replicate datasets are provided.",
         call. = FALSE)
  }

  if(method == "anova") {
    if(length(single.rep) == 0) {
      lcpm.aov <- aov(lcpm ~ samples, data = data.frame(samples = names(lcpm), lcpm = lcpm))
      aov.hsd <- TukeyHSD(lcpm.aov, conf.level=0.95)

      orig.comparisons <- rownames(aov.hsd$samples)
      # split.comparisons <- strsplit(orig.comparisons, "-")
      # new.comparisons <- sapply(
      #   X = split.comparisons, 
      #   FUN = function(x) {paste(x[2], "-", x[1], sep="")})

      comparisons <- data.frame(
        comparison = orig.comparisons,
        diff = aov.hsd$samples[,1],
        lower = aov.hsd$samples[,3],
        upper = aov.hsd$samples[,2],
        p.val = aov.hsd$samples[,4],
        p.adj = aov.hsd$samples[,4]
      )

      rownames(comparisons) <- orig.comparisons
    } else{
      stop("Cannot perform ANOVA if all data are not replicated.", call. = FALSE)
    }
  } else {
    comparisons <- list()

    if(length(multi.rep) > 1){
      for(i in seq(1, (length(multi.rep) - 1))){
        for(j in seq((i + 1), length(multi.rep))){
          if(multi.rep[i] != multi.rep[j]){
            s1 <- lcpm[names(lcpm) == multi.rep[j]]
            s2 <- lcpm[names(lcpm) == multi.rep[i]]
            contrast.name <- paste(multi.rep[j], "-", multi.rep[i], sep = "")
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

    if(length(single.rep) != 0){
      for(i in seq(1, length(multi.rep))){
        for(j in seq(1, length(single.rep))){
          if(multi.rep[i] != single.rep[j]){
            s1 <- lcpm[names(lcpm) == multi.rep[j]]
            s2 <- lcpm[names(lcpm) == single.rep[i]]
            contrast.name <- paste(multi.rep[j], "-", single.rep[i], sep = "")

            m1 <- mean(s1)
            n1 <- length(s1)
            m2 <- mean(s2)
            n2 <- length(s2)

            se <- sqrt((1 / n1 + 1 / n2) * ((n1 - 1) * sd(s1)^2 + (n2 - 1) * 0^2) / (n1 + n2 - 2))
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

      if(length(single.rep) > 1){
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
