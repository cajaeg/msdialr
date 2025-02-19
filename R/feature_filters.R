#' Subtract blank levels
#'
#' Summarize blank levels from individual blank samples using function 'fun' and
#' subtract summarized values from sample intensities. Possible negative values
#' resulting from subtraction are replaced by zero.
#' @param x MS-DIAL alignment results table
#' @param blank_idx index of blank samples; \code{NULL} to use \code{attr(x,
#'   "msdial_sam")$file_type == "Blank"}
#' @param fun function to use for blank level calculation, default: "max"
#'
#' @returns an object of the same class as 'x'
#' @export
subtractBlankLevels <- function(x, 
                                blank_idx = NULL, 
                                fun = c("max", "median", "average")[1]) {
  if (is.null(blank_idx)) {
    blank_idx <- attr(x, "msdial_sam")$file_type == "Blank"
  }
  int_mat <- data.matrix(getIntensityMatrix(x))
  int_blank <- pmax(0, apply(int_mat[, blank_idx], 1, fun), na.rm = TRUE)
  int_mat_new <- matrix(
    pmax(0, int_mat - int_blank, na.rm = TRUE),
    nrow = nrow(int_mat),
    ncol = ncol(int_mat)
  )
  x[, getIntensityColumns(x)] <- int_mat_new
  x
}


#' Remove features occurring only in few samples
#'
#' @param x MS-DIAL alignment results table
#' @param sample_group numeric or character vector defining sample grouping;
#'   \code{NULL} to use \code{attr(x, "msdial_sam")$class}; "" to apply filter
#'   without sample grouping (equivalent to \code{rep(1, length(samples))}).
#' @param min_int minimum intensity that a feature must have, default: 1e3
#' @param min_frac minimum fraction of samples (of at least one group if
#'   sample_group is defined accordingly) that contain the feature at an intensity >= min_int, default: 0.5
#'
#' @returns an object of the same class as 'x'
#' @export
removeSparseFeatures <- function(x,
                                 sample_group = NULL,
                                 min_int = 1e3,
                                 min_frac = 0.5) {
  stopifnot(inherits(x, "data.frame"))
  if (is.null(sample_group) && checkSam(x))
    sample_group <- attr(x, "msdial_sam")$class
  if(is.null(sample_group) || all(sample_group == ""))
    sample_group <- rep(1, length(getIntensityColumns(x)))
  int_mat <- getIntensityMatrix(x)
  flt <- apply(int_mat, 1, tapply, sample_group, function(x)
    sum(x >= min_int) / length(x) >= min_frac
  )
  if(!is.null(dim(flt))) 
    flt <- apply(flt, 2, any)
  x[flt, ]
}
