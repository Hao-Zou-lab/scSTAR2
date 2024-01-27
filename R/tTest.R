#' pearson correlation significant test
#'
#' @param r Correlation coefficient.
#' @param n Sample size.
#'
#' @return median_P,The tTest function is used to calculate the significance test of the Pearson correlation.


tTest <- function(r, n) {
  p <- pbeta(1 - r^2, n / 2 - 1, 1 / 2)
}

