#' Global Overlap Comparisons Effect Sizes
#'
#' Returns ROC AUC, Cohen's h, and Jensen-Shannon divergence values for feature overlap between conditions.
#' The reported AUC values are directionless -- always between 0.5 and 1.
#'
#' @param xint.obj An xIntOverlap object containing feature overlap information
#' @return A data frame with effect size statistics (ROC AUC and Jensen-Shannon divergence) for each comparison
#'
#' @importFrom pROC roc auc ci.auc
#' @importFrom stats qnorm
#'
#' @export
#'
overlap_effect <- function(xint.obj) {
  if(!validObject(xint.obj)) {
    stop("xint.obj is not a valid xIntObject.", call. = FALSE)
  }

  col.dat <- colData(xint.obj)
  conditions <- levels(col.dat$condition)

  cond.pairs <- combn(conditions, 2, simplify = FALSE)

  results <- lapply(
    X = seq_along(cond.pairs),
    FUN = function(i) {
      cond.pair <- cond.pairs[[i]]
      cond1 <- cond.pair[2]
      cond2 <- cond.pair[1]

      samples1 <- col.dat[col.dat$condition == cond1, ]$sample
      samples2 <- col.dat[col.dat$condition == cond2, ]$sample

      o.counts1 <- col.dat[samples1, "overlapping.sites"]
      t.counts1 <- col.dat[samples1, "total.sites"]
      o.counts2 <- col.dat[samples2, "overlapping.sites"]
      t.counts2 <- col.dat[samples2, "total.sites"]

      ts.1 <- sum(t.counts1)
      os.1 <- sum(o.counts1)
      score.1 <- c(rep(1, os.1), rep(0, ts.1 - os.1))
      ts.2 <- sum(t.counts2)
      os.2 <- sum(o.counts2)
      score.2 <- c(rep(1, os.2), rep(0, ts.2 - os.2))

      labels <- c(rep(0, ts.1), rep(1, ts.2))
      scores <- c(score.1, score.2)

      roc <- roc(response = labels, predictor = scores, quiet = TRUE)
      auc <- auc(roc)

      # Flip AUC value if < 0.5
      if (auc < 0.5) {
        auc <- 1 - auc
      }

      # prop1 <- sum(o.counts1) / sum(t.counts1)
      # prop2 <- sum(o.counts2) / sum(t.counts2)
      # phi1 <- 2 * asin(sqrt(prop1))
      # phi2 <- 2 * asin(sqrt(prop2))
      # h <- abs(phi1 - phi2)

      rate.scale <- 10

      if (length(o.counts1) == 1 || length(o.counts2) == 1) {
        # If no replication, fall back to Poisson
        lambda1 <- sum(o.counts1) / sum(t.counts1) * rate.scale
        lambda2 <- sum(o.counts2) / sum(t.counts2) * rate.scale
      } else {
        # If replication, use Quasi-Poisson
        qp1 <- glm(o.counts1 ~ t.counts1, family = quasipoisson())
        qp2 <- glm(o.counts2 ~ t.counts2, family = quasipoisson())

        lambda1 <- sum(predict(qp1, type = "response")) / sum(t.counts1) * rate.scale
        lambda2 <- sum(predict(qp2, type = "response")) / sum(t.counts2) * rate.scale
      }

      kl.1to2 <- lambda1 - lambda2 + lambda2 * log(lambda2 / lambda1)
      kl.2to1 <- lambda2 - lambda1 + lambda1 * log(lambda1 / lambda2)

      js <- 0.5 * (kl.1to2 + kl.2to1)

      out <- data.frame(
        comparison = paste0(cond1, "-", cond2),
        roc.auc = auc,
        # cohens.h = h,
        js.divergence = js
      )

      return(out)
    }
  )

  return(do.call(rbind, results))
}