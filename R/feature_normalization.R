#' Divide intensity values by summed or median intensity
#'
#' Normalization is performed by (1) calculating each sample's summed or median
#' intensity; (2) summarizing these values as group medians; (3) dividing
#' sample-wise summary values by group medians; (4) dividing peak intensities by
#' these normalized values. This way, normalized intensity remain similar
#' in magnitude as the original values.
#' @param x MS-DIAL alignment results table
#' @param sample_group numeric or character vector defining sample grouping;
#'   \code{NULL} to use \code{attr(x, "msdial_sam")$class}; "" to apply filter without sample
#'   grouping (equivalent to \code{rep(1, length(samples))})
#' @param fun summary function to use, i.e. sum or median
#'
#' @returns an object of the same class as 'x'
#' @export
normalizeIntensities <- function(x,
                                 sample_group = NULL,
                                 fun = c("sum", "median")[1]) {
  if(is.null(sample_group) && checkSam(x))
    sample_group <- attr(x, "msdial_sam")$class
  if(is.null(sample_group) || all(sample_group == ""))
    sample_group <- rep(1, length(getIntensityColumns(x)))
  int_mat <- getIntensityMatrix(x)
  stat <- apply(int_mat, 2, fun)
  stat_med <- stats::ave(stat, sample_group, FUN = stats::median)
  stat <- stat / stat_med
  int_mat_new <- matrix(pmax(0, t(t(int_mat) / stat), na.rm = TRUE),
                        nrow = nrow(int_mat), ncol = ncol(int_mat))
  x[, getIntensityColumns(x)] <- int_mat_new
  x
}
