#' Determine outliers
#' 
#' Determine outliers of a numeric vector, based on the first and third quartile minus/plus 1.5 times the
#' inter-quartile range.
#' @param x numeric vector containing values test for outliers
#' @export
#' @return a boolean vector, containing TRUE for outlier elements and FALSE for all other elements.
#' @usage \code{is_outlier(x)}
is_outlier <- function(x) {
  (x < quantile(x, 0.25) - 1.5 * IQR(x)) | (x > quantile(x, 0.75) + 1.5 * IQR(x))
}