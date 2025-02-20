#' Subtract blank levels
#'
#' Summarize blank levels from individual blank samples using function 'fun' and
#' subtract summarized values from sample intensities. Possible negative values
#' resulting from subtraction are replaced by zero.
#' @param x MS-DIAL alignment results table
#' @param blank_idx index of blank samples; \code{NULL} to use \code{attr(x,
#'   "msdial_sam")$file_type == "Blank"}
#' @param fun function to use for blank level calculation, default: "max"
#' @param ... further arguments passed to 'fun'
#'
#' @returns an object of the same class as 'x'
#' @export
subtractBlankLevels <- function(x, 
                                blank_idx = NULL, 
                                fun = c("max", "median", "average")[1],
                                ...) {
  if (is.null(blank_idx)) {
    blank_idx <- attr(x, "msdial_sam")$file_type == "Blank"
  }
  int_mat <- getIntensityMatrix(x, as.matrix = TRUE)
  int_blank <- pmax(0, apply(int_mat[, blank_idx], 1, fun, ...))
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
#' @param verbose display messages
#'
#' @returns an object of the same class as 'x'
#' @export
removeSparseFeatures <- function(x,
                                 sample_group = NULL,
                                 min_int = 1e3,
                                 min_frac = 0.5,
                                 verbose = NULL) {
  stopifnot(inherits(x, "data.frame"))
  if(is.null(verbose))
    verbose <- interactive()
  msg <- function(x, verbose) 
    if(verbose)
      message(x)
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
  x_out <- x[flt, ]
  msg(sprintf("Number of peaks reduced from %d to %d", 
              nrow(x), 
              nrow(x_out)),
      verbose)
  x_out
}


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
#' @param ... further arguments passed to 'fun'
#'
#' @returns an object of the same class as 'x'
#' @export
normalizeIntensities <- function(x,
                                 sample_group = NULL,
                                 fun = c("sum", "median")[1],
                                 ...) {
  if(is.null(sample_group) && checkSam(x))
    sample_group <- attr(x, "msdial_sam")$class
  if(is.null(sample_group) || all(sample_group == ""))
    sample_group <- rep(1, length(getIntensityColumns(x)))
  int_mat <- getIntensityMatrix(x, as.matrix = TRUE)
  stat <- apply(int_mat, 2, fun, ...)
  stat_med <- stats::ave(stat, sample_group, FUN = stats::median)
  stat <- stat / stat_med
  int_mat_new <- matrix(pmax(0, t(t(int_mat) / stat), na.rm = TRUE),
                        nrow = nrow(int_mat), ncol = ncol(int_mat))
  x[, getIntensityColumns(x)] <- int_mat_new
  x
}
