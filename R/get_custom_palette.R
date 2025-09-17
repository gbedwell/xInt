#' Custom Color Palette
#'
#' Creates a custom color palette for use in plotting functions. The palette is designed to provide visually distinct colors for up to 12 categories, with additional colors generated dynamically if more categories are needed.
#'
#' @param n Integer. The number of colors to generate. If \code{n} is less than or equal to 12, predefined base colors are used. If \code{n} is greater than 12, additional colors are interpolated using \code{colorRampPalette}.
#'
#' @return A character vector of color codes. The length of the vector matches the value of \code{n}.
#'
#' @examples
#' # Generate a palette for 5 categories
#' get_custom_palette(5)
#'
#' # Generate a palette for 15 categories
#' get_custom_palette(15)
#'
#' @export
get_custom_palette <- function(n) {
  base.colors <- c("darkred", "darkblue", "darkgreen", "darkorange", 
                   "#3B2042", "darkgoldenrod", "darkcyan", "darkmagenta",
                   "chocolate4", "darkslategray", "midnightblue", "saddlebrown")
  
  if (n <= length(base.colors)) {
    return(base.colors[1:n])
  } else {
    return(colorRampPalette(base.colors)(n))
  }
}