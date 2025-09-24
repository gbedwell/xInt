#' Plot Density
#'
#' Plot per-condition average feature density.
#'
#' @param x The output of feature_density() with options average = TRUE and conditions are given.
#' @param print.plot Boolean. Whether or not to print the plot(s). Defaults to FALSE.
#' @param plot.title The title of the plot. When NULL, no title is added to the plot.
#' 
#' @return A ggplot2 object.
#'
#' @import ggplot2
#'
#' @export
#'
plot_density <- function(x, stats.df = NULL, condition.levels = NULL, print.plot = FALSE, plot.title = NULL) {

  if(!all(c("sample", "condition", "avg.density") %in% colnames(x))) {
    stop("x must match the output of feature_density().",
         call. = FALSE)
  }

  if(!is.null(stats.df)) {
    if(!all(c("comparison", "p.adj") %in% colnames(stats.df))) {
      stop("stats.df must have columns 'comparison' and 'p.adj'.",
          call. = FALSE)
    }
  }

  print("This function is not yet completed!")
}



