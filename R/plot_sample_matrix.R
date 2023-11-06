#' Pairwise count matrix plot
#'
#' Plots summary information about each of the datasets used to build the xInt model.
#'
#'@param xint.obj The output of make_xInt_dataset().
#'@param labels The desired sample labels. Length should equal the number of columns in the count matrix.
#'
#'
#'@export
#'
plot_sample_matrix <- function (xint.obj, labels = NULL)
{
  frac.data <- SummarizedExperiment::colData(xint.obj)$fraction.overlap
  count.data <- SummarizedExperiment::assays(xint.obj)$counts
  df <- data.frame(rbind(total = SummarizedExperiment::colData(xint.obj)$total.sites,
                         fraction = SummarizedExperiment::colData(xint.obj)$fraction.overlap,
                         SummarizedExperiment::assays(xint.obj)$counts))
  if (is.null(labels)) {
    labels = colnames(xint.obj)
  }
  if (length(labels) != ncol(count.data)) {
    stop("The number of labels must equal the number of columns in the count matrix.",
         call. = FALSE)
  }
  cor_func <- function(data, mapping, method, ...) {
    x <- GGally::eval_data_col(data[rownames(data) != "fraction" &
                                      rownames(data) != "total", ], mapping$x)
    y <- GGally::eval_data_col(data[rownames(data) != "fraction" &
                                      rownames(data) != "total", ], mapping$y)

    corr <- cor(x, y, method = method, use = "everything")

    calc.jac <- function(x, y) {
      both <- apply(X = matrix(c(x, y), ncol = 2), MARGIN = 1,
                    FUN = function(x) {
                      all.vec <- all(x > 0)
                      sum(all.vec)
                    })
      one <- apply(X = matrix(c(x, y), ncol = 2), MARGIN = 1,
                   FUN = function(x) {
                     one.vec <- any(x == 0) & !all(x == 0)
                     sum(one.vec)
                   })
      none <- apply(X = matrix(c(x, y), ncol = 2), MARGIN = 1,
                    FUN = function(x) {
                      none.vec <- all(x == 0)
                      sum(none.vec)
                    })
      (sum(both) + sum(none))/(sum(both) + sum(one) + sum(none))
      # sum(both)/(sum(both) + sum(one))
    }
    jac <- calc.jac(x, y)
    GGally::ggally_text(label = paste0("r = ", as.character(sprintf("%.3f",
                                                                    round(x = corr, digits = 3))),
                                       "\n",
                                       "J = ", as.character(sprintf("%.3f",
                                                                    round(x = jac, digits = 3)))),
                        mapping = ggplot2::aes(),
                        xP = 0.5, yP = 0.5, color = "darkred", ...) +
      ggplot2::theme_bw()
  }
  point_alpha <- function(data, mapping, ...) {
    ggplot2::ggplot(data = log(data[rownames(data) != "fraction" &
                                      rownames(data) != "total", ] + 1), mapping = mapping) +
      ggplot2::geom_point(..., shape = 21, size = 2, alpha = 0.15,
                          fill = "black") + ggplot2::theme_bw()
  }
  perc_int <- function(data, mapping, ...) {
    n.sites <- data[rownames(data) == "total", ]
    frac.in <- data[rownames(data) == "fraction", ]
    mapped <- GGally::eval_data_col(frac.in, mapping$x)
    mapped1 <- GGally::eval_data_col(n.sites, mapping$x)

    GGally::ggally_text(label = paste0(paste0("n = ",
                                              format(mapped1, big.mark = ",", scientific = FALSE)),
                                       "\n",
                                       paste0("(",
                                              as.character(sprintf("%.2f",
                                                                   round(x = mapped * 100, digits = 2))), "%)")),
                        mapping = ggplot2::aes(),
                        xP = 0.5, yP = 0.5, color = "darkblue", ...) + ggplot2::theme_bw() +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank())
  }
  GGally::ggpairs(df, upper = list(continuous = GGally::wrap(cor_func,
                                                             method = "pearson")), lower = list(continuous = GGally::wrap(point_alpha)),
                  diag = list(continuous = GGally::wrap(perc_int, size = 3.5)), columnLabels = labels,
                  xlab = "log(Counts + 1)", ylab = "log(Counts + 1)", progress = FALSE) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "gray85",
                                                            colour = "black", size = 0.5), axis.title = ggplot2::element_text(size = 14))
}
