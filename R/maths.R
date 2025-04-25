# project: omx
# XQ - LUMC
# april 25 - 2025
# math notation for plotting


#' Format Scientific Notation for Plot Labels
#'
#' Converts a numeric value into a formatted scientific notation expression suitable for plot labels.
#'
#' @param x Numeric value to be formatted.
#'
#' @return A parsed expression representing the value in scientific notation (e.g., 1.23 * 10^4).
#'
#' @examples
#' \dontrun{
#'   format_scientific_expr(0.000123)  # Returns expression for 1.23 * 10^-4
#' }
format_scientific_expr <- function(x) {
  sci <- format(x, scientific = TRUE)
  coef <- as.numeric(sub("e.*", "", sci))
  exp <- as.numeric(sub(".*e", "", sci))
  parse(text = sprintf("%g %%*%% 10^%d", coef, exp))
}