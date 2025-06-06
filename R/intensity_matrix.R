#' Subtract blank levels
#'
#' Summarize blank levels from individual blank samples using function 'fun' and
#' subtract summarized values from sample intensities. Possible negative values
#' resulting from subtraction are replaced by zero.
#' @param x MS-DIAL alignment results table
#' @param blank_idx numeric or logical index of blank samples; if \code{NULL},
#'   \code{getSampleList(x)$file_type == "Blank"} is used. If this evaluates to
#'   NULL, 'x' is returned unchanged.
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
    blank_idx <- getSampleList(x)$file_type == "Blank"
  }
  if(!is.null(blank_idx)) {
    if((is.logical(blank_idx) && sum(blank_idx) > 0) ||
       (is.numeric(blank_idx) && length(blank_idx) > 0)) {
      int_mat <- getIntensityMatrix(x, as.matrix = TRUE)
      int_blank <- pmax(0, apply(int_mat[, blank_idx], 1, fun, ...))
      int_mat_mod <- matrix(
        pmax(0, int_mat - int_blank, na.rm = TRUE),
        nrow = nrow(int_mat),
        ncol = ncol(int_mat)
      )
      x[, getIntensityColumns(x)] <- int_mat_mod
      x
    } else {
      warning("blank index does not indicate any blank sample")
      x
    }
  } else {
    warning("blank index is NULL")
    x
  }
}


#' Flag/filter features occurring only in few samples
#'
#' @param x MS-DIAL alignment results table
#' @param sample_group numeric or character vector defining sample grouping;
#'   \code{NULL} to use \code{attr(x, "msdial_sam")$class}; "" to apply filter
#'   without sample grouping (equivalent to \code{rep(1, length(samples))}).
#' @param min_int minimum intensity that a feature must have, default: 1e3
#' @param min_frac minimum fraction of samples (of at least one group if
#'   sample_group is defined accordingly) that contain the feature at an
#'   intensity >= min_int, default: 0.5
#' @param flag_column name of flag column, default: "passes_sparse_filter"
#' @param filter filter results according to flag, default: FALSE
#' @param verbose display messages, defaults to the output of
#'   \code{interactive()}
#'
#' @returns an object of the same class as 'x'
#' @export
filterSparseFeatures <- function(x,
                                 sample_group = NULL,
                                 min_int = 1e3,
                                 min_frac = 0.5,
                                 flag_column = "passes_sparse_filter",
                                 filter = FALSE,
                                 verbose = NULL) {
  stopifnot(inherits(x, "data.frame"))
  class_column <- "class"
  if(is.null(verbose))
    verbose <- interactive()
  msg <- function(x, verbose) 
    if(verbose)
      message(x)
  if (is.null(sample_group))
    sample_group <- getSampleList(x)[[class_column]]
  if(is.null(sample_group) || all(sample_group == ""))
    sample_group <- rep(1, length(getIntensityColumns(x)))
  int_mat <- getIntensityMatrix(x)
  stopifnot(ncol(int_mat) == length(sample_group))
  flag <- apply(int_mat, 1, tapply, sample_group, function(x)
    sum(x >= min_int) / length(x) >= min_frac
  )
  if(!is.null(dim(flag))) 
    flag <- apply(flag, 2, any)
  x[[ flag_column ]] <- flag
  x_out <- if(filter) {
    x[flag, ]
  } else {
    x
  }
  msg(sprintf("Found %d of %d peaks passing sparse feature filter", 
              sum(flag), 
              nrow(x)),
      verbose)
  x_out |> 
    relocateIntensityColumns()
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
  class_column = "class"
  if(is.null(sample_group))
    sample_group <- getSampleList(x)[[class_column]]
  if(is.null(sample_group) || all(sample_group == ""))
    sample_group <- rep(1, length(getIntensityColumns(x)))
  int_mat <- getIntensityMatrix(x, as.matrix = TRUE)
  stat <- apply(int_mat, 2, fun, ...)
  stat_med <- stats::ave(stat, sample_group, FUN = stats::median)
  stat <- stat / stat_med
  int_mat_mod <- matrix(pmax(0, t(t(int_mat) / stat), na.rm = TRUE),
                        nrow = nrow(int_mat), ncol = ncol(int_mat))
  x[, getIntensityColumns(x)] <- int_mat_mod
  x
}


#' Filter feature table according to applied flag/rank functions
#'
#' @param x MS-DIAL alignment results table
#' @param which index selecting filter columns from \code{c("passes_groupsize_filter",
#'   "rank_MS2", "passes_sparse_filter")}, default: all (\code{c(1:3)})
#'
#' @returns filtered alignment table
#' @noRd
#' @export
applyFilters <- function(x, which = 1:3) {
  col_names <- c("passes_groupsize_filter",
                 "rank_MS2",
                 "passes_sparse_filter")[which]
  stopifnot(all(col_names %in% colnames(x)))
  flt <- data.frame(x[, col_names])
  if("rank_MS2" %in% col_names) {
    flt[, "rank_MS2"] <- is.finite(flt[, "rank_MS2"]) & flt[, "rank_MS2"] == 1
  }
  flt <- apply(flt, 1, all)
  x[flt, ]
}