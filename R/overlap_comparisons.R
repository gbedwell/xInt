#' GLM Feature Overlap Comparisons
#'
#' Perform per-feature or feature-set integration targeting comparisons using (quasi-)Poisson GLMs.
#' Supports Poisson and Quasi-Poisson models with flexible formula specification.
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
#' @param formula A formula object or character string specifying the model. If NULL (default),
#'   a formula will be constructed based on the condition variable and any provided batch variables.
#' @param batch A character vector of column names in the colData of the xIntObject containing relevant
#'   batch information. Defaults to NULL.
#' @param contrasts Character vector of contrasts to test, in the form "A-B", denoting pairs and directionality. 
#' If NULL (default), all pairwise comparisons between conditions will be created.
#' @param type Character string specifying the type of comparisons: per-feature ('local') or feature-set ('global').
#' Default 'local'.
#' @param intercept Boolean. Whether or not to include an intercept term in the model. 
#' Including an intercept means coefficient effects and significance values are defined relative to a reference
#' (i.e., the first term in conditions).
#' Default TRUE.
#' @param model.family The model family to use in the GLM. Can be NULL (default), 'poisson', or 'quasipoisson'.
#' When NULL, model.family is set to 'poisson' when type = 'local' and 'quasipoisson' when type = 'global'.
#' 
#' @return A list containing the results of GLM fits for each feature, including:
#'   \item{kept.feats}{Vector of features that passed filtering}
#'   \item{perc.change}{Matrix of percentage changes for pairwise comparisons}
#'   \item{log2.fc}{Matrix of log2 fold changes for pairwise comparisons}
#'   \item{coefficients}{Matrix of coefficient estimates}
#'   \item{p.vals}{Matrix of p-values for each coefficient}
#'   \item{contrasts}{Matrix of log fold changes for pairwise comparisons}
#'   \item{contrast.p.vals}{Matrix of p-values for pairwise comparisons}
#'   \item{model.info}{Information about the model used}
#'   \item{contrast.mat}{The contrast matrix used for comparisons}
#'
#' @examples
#' data(xobj)
#' local_comparisons(xint.obj = xobj, plot = FALSE)
#'
#' @importFrom stats model.matrix glm poisson quasipoisson coef vcov pchisq pnorm
#' @import SummarizedExperiment
#' @import ggplot2
#'
#' @export
#'
overlap_comparisons <- function(xint.obj, min.count, min.total.count, 
                                remove.zeros = TRUE, plot = TRUE, formula = NULL, 
                                batch = NULL, contrasts = NULL, type = c("local", "global"), 
                                intercept = TRUE, model.family = NULL) {
  
  if(!validObject(xint.obj)) {
    stop("xint.obj is not a valid xIntObject.", call. = FALSE)
  }

  type <- match.arg(type)
  
  if(is.null(model.family)) {
    model.family <- ifelse(type == "local", "poisson", "quasipoisson")
  }
  
  if(!model.family %in% c("poisson", "quasipoisson")){
    stop("model.family must be one of 'poisson' or 'quasipoisson'.", call. = FALSE)
  }
  
  if(model.family == "poisson"){
    family.type = poisson()
  } else if(model.family == "quasipoisson"){
    family.type = quasipoisson()
  }
  
  dat <- colData(xint.obj)

  if(length(unique(dat$condition)) == nrow(dat)){
    stop("No replication detected. To compare unreplicated data, use single_overlap_comparisons().",
    call. = FALSE)
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
  
  # Create design matrix based on formula or default
  if(is.null(formula)) {
    if(!is.null(batch)) {
      if(intercept){
        formula.str <- paste("~ condition +", paste(batch, collapse = " + "))
      } else{
        formula.str <- paste("~ 0 + condition +", paste(batch, collapse = " + "))
      }
    } else {
      if(intercept){
        formula.str <- "~ condition"
      } else{
        formula.str <- "~ 0 + condition"
      }
    }
    formula <- as.formula(formula.str)
  } else if(is.character(formula)) {
    formula <- as.formula(formula)
  }
  
  design <- model.matrix(formula, data = dat)
  
  # Clean up column names
  colnames(design) <- gsub(pattern = "condition", replacement = "", x = colnames(design))
  
  if(!is.null(batch)) {
    for(entry in batch) {
      colnames(design) <- gsub(pattern = entry, replacement = "", x = colnames(design))
    }
  }
  
  # Filter low count features
  keep <- rowSums(cs >= min.count) >= min.total.count
  filtered.cs <- cs[keep, , drop = FALSE]

  lcpm <- calculate_cpm(filtered.cs, lib.size = dat$total.sites, log = TRUE, prior.count = 2)
  
  if(isTRUE(plot) && type == "local") {
    plotMDS(lcpm, labels = dat$condition)
    plot_meanvar(filtered.counts = filtered.cs, lib.sizes = dat$total.sites,
                dat = dat, design = design)  
  }
  
  # Get condition columns for contrasts
  # This is more complex with custom formulas, so we'll extract main effect terms
  main.terms <- attr(terms(formula), "term.labels")
  main.terms <- main.terms[!grepl(":", main.terms)]  # Remove interaction terms
  
  # Try to identify condition columns
  condition.cols <- colnames(design)
  if("(Intercept)" %in% condition.cols) {
    condition.cols <- condition.cols[condition.cols != "(Intercept)"]
    all.levels <- levels(dat$condition)
    ref.level <- all.levels[!all.levels %in% condition.cols][1]
    
    if(!is.null(ref.level)) {
      condition.cols <- c(ref.level, condition.cols)
    }
  }
  
  if(!is.null(batch)) {
    for(b in batch) {
      condition.cols <- condition.cols[!grepl(b, condition.cols)]
    }
  }
  
  # Create contrast matrix
  if(is.null(contrasts)) {
    # Create all pairwise contrasts
    contrast.names <- character()
    for(i in seq(1, (length(condition.cols) - 1))) {
      for(j in seq((i + 1), length(condition.cols))) {
        if(condition.cols[i] != condition.cols[j]) {
          contrast.name <- paste(condition.cols[j], "-", condition.cols[i], sep = "")
          contrast.names <- c(contrast.names, contrast.name)
        }
      }
    }
  } else {
    contrast.names <- contrasts
  }
  
  contrast.mat <- make_contrasts(
    design.cols = colnames(design),
    condition.cols = condition.cols,
    contrasts = contrast.names,
    has.intercept = "(Intercept)" %in% colnames(design)
  )
  
  # Calculate offsets
  log.offsets <- log(dat$total.sites)
  
  # Fit GLM for each feature
  n.features <- nrow(filtered.cs)
  n.samples <- ncol(filtered.cs)
  n.coef <- ncol(design)
  n.contrasts <- ncol(contrast.mat)
  
  coef.mat <- matrix(NA, nrow = n.features, ncol = n.coef)
  pval.mat <- matrix(NA, nrow = n.features, ncol = n.coef)
  contrast.mat.results <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  percent.change.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  log2.fc.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  contrast.pval.mat <- matrix(NA, nrow = n.features, ncol = n.contrasts)
  
  rownames(coef.mat) <- rownames(pval.mat) <- 
    rownames(contrast.mat.results) <- rownames(contrast.pval.mat) <- 
    rownames(percent.change.mat) <- rownames(log2.fc.mat) <- rownames(filtered.cs)
  
  colnames(coef.mat) <- colnames(pval.mat) <- colnames(design)
  colnames(contrast.mat.results) <- colnames(contrast.pval.mat) <- 
    colnames(percent.change.mat) <- colnames(log2.fc.mat) <- colnames(contrast.mat)
  
  for(i in 1:n.features) {
    # Fit GLM
    fit.glm <- glm.fit(
      x = design, 
      y = filtered.cs[i, ], 
      family = family.type, 
      offset = log.offsets
    )
    
    # Create glm object
    fit <- list(
      coefficients = fit.glm$coefficients,
      residuals = fit.glm$residuals,
      fitted.values = fit.glm$fitted.values,
      effects = fit.glm$effects,
      R = if(is.null(fit.glm$qr)) NULL else fit.glm$qr$qr,
      rank = fit.glm$rank,
      qr = fit.glm$qr,
      family = family.type,
      linear.predictors = fit.glm$linear.predictors,
      deviance = fit.glm$deviance,
      aic = fit.glm$aic,
      null.deviance = fit.glm$null.deviance,
      iter = fit.glm$iter,
      weights = fit.glm$weights,
      prior.weights = fit.glm$prior.weights,
      df.residual = fit.glm$df.residual,
      df.null = fit.glm$df.null,
      y = filtered.cs[i, ],
      x = design,
      model = data.frame(y = filtered.cs[i, ], as.data.frame(design)),
      converged = fit.glm$converged,
      boundary = fit.glm$boundary,
      call = call("glm"),
      formula = as.formula("y ~ ."),
      terms = terms(as.formula("y ~ ."), data = as.data.frame(design)),
      data = data.frame(y = filtered.cs[i, ], as.data.frame(design)),
      offset = log.offsets,
      control = list(),
      method = "glm.fit",
      contrasts = NULL,
      xlevels = list()
    )
    class(fit) <- c("glm", "lm")
    
    # Extract coefficients
    coef.mat[i, ] <- coef(fit)
    
    # Calculate p-values
    fit.summary <- summary(fit)
    pval.mat[i, ] <- fit.summary$coefficients[, 4]
    vcov.mat <- vcov(fit)
    
    # Apply contrasts using the contrast matrix
    for(j in 1:ncol(contrast.mat)) {
      # Get contrast vector
      contrast.vector <- contrast.mat[, j]
      
      # Calculate contrast estimate (similar to limma's contrasts.fit)
      contrast.estimate <- sum(coef(fit) * contrast.vector)
      contrast.mat.results[i, j] <- contrast.estimate
      
      # Calculate percent change
      percent.change.mat[i, j] <- (exp(contrast.estimate) - 1) * 100
      
      # Calculate fold change
      log2.fc.mat[i, j] <- contrast.estimate / log(2)
      
      # Calculate variance of contrast
      contrast.var <- t(contrast.vector) %*% vcov.mat %*% contrast.vector
      
      # Calculate p-value using Wald test
      wald.stat <- contrast.estimate^2 / contrast.var
      contrast.pval.mat[i, j] <- pchisq(wald.stat, df = 1, lower.tail = FALSE)
    }
  }
  
  # Return results
  results <- list(
    kept.feats = rownames(filtered.cs),
    log.cpm = lcpm,
    perc.change = percent.change.mat,
    log2.fc = log2.fc.mat,
    coefficients = coef.mat,
    p.values = pval.mat,
    contrasts = contrast.mat.results,
    contrast.p.vals = contrast.pval.mat,
    model.info = list(
      formula = formula,
      family = family.type,
      design = design
    ),
    contrast.mat = contrast.mat
  )

  if(isTRUE(plot) && type == "local") {
    plot_pvals(ps = c(results$contrast.p.vals))
  }
  
  return(results)
}

#' Create a contrast matrix for GLM comparisons
#'
#' @param levels Character vector of coefficient names from the design matrix
#' @param contrasts Character vector of contrast formulas (e.g., "A-B", "A-C")
#'
#' @return A contrast matrix
#'
make_contrasts <- function(design.cols, condition.cols, contrasts, has.intercept = FALSE) {
  n.design.cols <- length(design.cols)
  n.contrasts <- length(contrasts)
  
  contrast.mat <- matrix(0, nrow = n.design.cols, ncol = n.contrasts)
  rownames(contrast.mat) <- design.cols
  colnames(contrast.mat) <- contrasts
  
  for(i in 1:n.contrasts) {
    contrast.formula <- contrasts[i]
    parts <- strsplit(contrast.formula, "-")[[1]]
    
    if(length(parts) != 2) {
      stop("Contrast must be in the form 'A-B'")
    }
    
    cond1 <- parts[1]
    cond2 <- parts[2]
    
    if(has.intercept) {
      # Find the reference level (not in design.cols but in condition.cols)
      ref.level <- condition.cols[!condition.cols %in% design.cols]
      ref.level <- ref.level[!ref.level %in% c("(Intercept)")]
      
      # Automatically handle contrasts involving the reference level
      if(cond2 %in% ref.level) {
        # 1 * other - 1 * ref
        if(cond1 %in% design.cols) {
          contrast.mat[cond1, i] <- 1
          # Reference level is not in design.cols, so skip assigning to contrast.mat[cond1, i]
        } else {
          stop(paste("Condition", cond1, "not found in design matrix"), call. = FALSE)
        }
      } else if(cond1 %in% ref.level) {
        if(cond2 %in% design.cols) {
          contrast.mat[cond2, i] <- 1
          # Reference level is not in design.cols, so skip assigning to contrast.mat[cond1, i]
        } else {
          stop(paste("Condition", cond2, "not found in design matrix"), call. = FALSE)
        }
      } else {
        # 1 * cond1 - 1 * cond2
        if(cond1 %in% design.cols) {
          contrast.mat[cond1, i] <- 1
        } else {
          stop(paste("Condition", cond1, "not found in design matrix"), call. = FALSE)
        }
        
        if(cond2 %in% design.cols) {
          contrast.mat[cond2, i] <- -1
        } else {
          stop(paste("Condition", cond2, "not found in design matrix"), call. = FALSE)
        }
      }
    } else {
      if (!(cond1 %in% design.cols) || !(cond2 %in% design.cols)) {
        stop("Contrast conditions must be in the design matrix")
      }
      
      contrast.mat[cond1, i] <- 1
      contrast.mat[cond2, i] <- -1
    }
  }
  
  return(contrast.mat)
}

plot_meanvar <- function(filtered.counts, lib.sizes, dat = NULL, design = NULL){  
  # Calculate offset-adjusted means and variances
  n.features <- nrow(filtered.counts)
  means <- numeric(n.features)
  vars <- numeric(n.features)

  # Use design matrix if provided, otherwise fall back to intercept-only
  if(!is.null(design) && !is.null(dat)) {
    # Use the actual experimental design
    for (i in 1:n.features) {
      # Fit GLM with full design matrix (same as main analysis)
      fit <- glm(filtered.counts[i,] ~ ., data = as.data.frame(design[, -1, drop = FALSE]), 
                 family = poisson(), offset = log(lib.sizes))
      fitted.values <- predict(fit, type = "response")
      pearson.residuals <- (filtered.counts[i,] - fitted.values) / sqrt(fitted.values)    
      means[i] <- mean(fitted.values)
      vars[i] <- means[i] * mean(pearson.residuals^2)
    }
  } else if(!is.null(dat) && "condition" %in% colnames(dat)) {
    # Use condition variable if available
    for (i in 1:n.features) {
      fit <- glm(filtered.counts[i,] ~ condition, data = dat, 
                 family = poisson(), offset = log(lib.sizes))
      fitted.values <- predict(fit, type = "response")
      pearson.residuals <- (filtered.counts[i,] - fitted.values) / sqrt(fitted.values)    
      means[i] <- mean(fitted.values)
      vars[i] <- means[i] * mean(pearson.residuals^2)
    }
  } else {
    # Fall back to original intercept-only approach
    for (i in 1:n.features) {
      fit <- glm(filtered.counts[i,] ~ 1, family = poisson(), offset = log(lib.sizes))    
      fitted.values <- predict(fit, type = "response")
      pearson.residuals <- (filtered.counts[i,] - fitted.values) / sqrt(fitted.values)    
      means[i] <- mean(fitted.values)
      vars[i] <- means[i] * mean(pearson.residuals^2)
    }
  }

  # Calculate dispersion values (variance/mean ratio)
  dispersions <- vars / means
  dispersion.ratio <- mean(dispersions, na.rm = TRUE)
  median.dispersion <- median(dispersions, na.rm = TRUE)
  
  plot.data <- data.frame(
    mean = means,
    variance = vars,
    log.mean = log10(means + 0.5),
    log.var = log10(vars + 0.5)
  )
  
  n.bins <- 50
  
  sorted.data <- plot.data[order(plot.data$log.mean),]
  
  # Create density-based bins
  points.per.bin <- floor(nrow(sorted.data) / n.bins)
  bin.indices <- rep(1:n.bins, each = points.per.bin)
  
  # Handle remaining points
  remaining <- nrow(sorted.data) - length(bin.indices)
  if(remaining > 0) {
    bin.indices <- c(bin.indices, rep(n.bins, remaining))
  }
  
  sorted.data$bin <- bin.indices
  
  # Calculate mean and variance for each density-based bin
  binned.data <- stats::aggregate(
    cbind(log.mean, log.var) ~ bin, 
    data = sorted.data, 
    FUN = mean
  )

  x.range <- range(plot.data$log.mean, na.rm = TRUE)
  x.seq <- seq(min(x.range), max(x.range), length.out = 100)
  poisson.data <- data.frame(
    x = x.seq,
    y = x.seq,
    type = "Poisson"
  )
  
  # Add type to binned data
  binned.data$type <- "Binned Variance"
  
  # Update subtitle to indicate which model was used
  model.type <- if(!is.null(design)) {
    "Design-adjusted Dispersion"
  } else if(!is.null(dat) && "condition" %in% colnames(dat)) {
    "Condition-adjusted Dispersion"
  } else {
    "Intercept-only Dispersion"
  }
  
  p <- ggplot() +
    geom_point(
      data = plot.data,
      aes(x = log.mean, y = log.var),
      color = "gray50", 
      alpha = 0.5,
      size = 2
    ) +
    geom_line(
      data = poisson.data,
      aes(x = x, y = y, color = type),
      linewidth = 1
    ) +
    geom_point(
      data = binned.data,
      aes(x = log.mean, y = log.var, color = type),
      size = 3
    ) +
    geom_line(
      data = binned.data,
      aes(x = log.mean, y = log.var, color = type),
      linewidth = 1
    ) +
    # Set colors for legend
    scale_color_manual(
      values = c("Poisson" = "darkred", "Binned Variance" = "darkblue"),
      name = NULL
    ) +
    # Add labels and title
    labs(
      title = "Mean-Variance Relationship",
      subtitle = paste("Average Dispersion:", round(dispersion.ratio, 2)),
      x = "log10(Mean + 0.5)",
      y = "log10(Variance + 0.5)"
    ) +
    # Use a clean theme
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      legend.justification = "left",
      legend.box.margin = margin(-10,-10,-10,-10),
      axis.text = element_text(size=12),
      axis.title = element_text(size=14)
    )
  
  # Print the plot
  print(p)
  
  # Return the plot object invisibly
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