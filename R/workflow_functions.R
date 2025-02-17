#' Group features based on identical neutral mass
#'
#' @param x MS-DIAL alignment results table
#' @param mz_column name of m/z column
#' @param adduct_column name of adduct column
#' @param rt_column name of retention time column
#' @param delta_mz absolute mass window (Da)
#' @param delta_rt absolute retention time window (min)
#' @param verbose display progress messages
#'
#' @returns 'x' with additional columns 'feat_group' and 'group_size'
#' @export
#' 
#' @examples
#' fp <- system.file("extdata/MSDIAL_Alignment_result_LC-MS.txt", package = "msdialr")
#' aligned <- loadAlignmentResults(fp)
#' length(unique(aligned$alignment)) # 230
#' 
#' aligned1 <- aligned |> 
#'   assignAdductGroups()
#' length(unique(aligned1$feat_group)) # 174
#' plotAdductSpace(aligned1)
#' 
#' aligned2 <- aligned1 |> 
#'   filterAdductGroupsBySize(min_group_size = 2)
#' length(unique(aligned2$feat_group)) # 36
#' plotAdductSpace(aligned2)
assignAdductGroups <- function(x, 
                               mz_column = "average_mz", 
                               adduct_column = "adduct_type",
                               rt_column = "average_rt_min",
                               delta_mz = 0.01,
                               delta_rt = 0.1, 
                               verbose = NULL) {
  if(is.null(verbose))
    verbose <- nrow(x) > 5e3
  msg <- function(x, verbose) 
    if(verbose)
      message(x)
  length_unique <- function(x) 
    length(unique(x))
  mz_digits <- 3
  rt_digits <- 2
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c(mz_column, adduct_column, rt_column) %in% colnames(x)))
  mz <- x[[ mz_column ]]
  adduct_type <- x[[ adduct_column ]]
  rt <- x[[ rt_column ]]
  msg("Calculating neutral masses ...", verbose)
  nm <- mz2nm(mz, adduct_type)
  msg("Grouping neutral masses ...", verbose)
  nm_group <- hcgroup(nm, h = delta_mz)
  msg("Grouping retention times ...", verbose)
  rt_group <- hcgroup(rt, h = delta_rt)
  feat_group <- paste(nm_group, rt_group, sep = "_")
  nm_group_lbl <- round(stats::ave(nm, feat_group, FUN = stats::median), mz_digits)
  rt_group_lbl <- round(stats::ave(rt, feat_group, FUN = stats::median), rt_digits)
  feat_group_lbl <- sprintf("%.3f_%.2f", nm_group_lbl, rt_group_lbl)
  group_size <- as.integer(stats::ave(adduct_type, feat_group, FUN = length_unique))
  x$feat_group <- feat_group_lbl
  x$group_size <- group_size
  x |> 
    relocateIntensityColumns() |> 
    updateIntensityColumnIndex()
}


#' Plot adduct groups
#' @rdname assignAdductGroups
#'
#' @param x MS-DIAL alignment results table
#' @param rt_column name of retention time column
#' @param mz_column name of m/z column
#' @param adduct_column name of adduct column
#'
#' @returns NULL
#' @export
plotAdductGroups <- function(x, 
                            rt_column = "average_rt_min",
                            mz_column = "average_mz",
                            adduct_column = "adduct_type") {
  rt <- x[[ rt_column ]]
  mz <- x[[ mz_column ]]
  fct <- factor(x[[ adduct_column ]])
  tmp <- sort(table(fct), decreasing = TRUE)
  fct <- factor(fct, levels = levels(fct)[match(names(tmp), levels(fct))])
  bg_col <- match(fct, levels(fct))
  plot(rt, mz, type = "n", 
       xlab = "Retention Time (min)", ylab = "m/z")
  points(rt, mz, 
         col = NULL, pch = 21, bg = bg_col)
  legend("topleft", legend = sprintf("%s (%d)", levels(fct), as.numeric(tmp)), 
         pch = 21, col = NA, pt.bg = seq_along(levels(fct)))
}


#' Remove features with less than n adducts
#'
#' @param x MS-DIAL alignment results table
#' @param group_size_column name of group size column, default: "group_size"
#' @param min_group_size minimum number of adducts
#'
#' @returns an object of the same class as 'x'
#' @export
filterAdductGroupsBySize <- function(x, 
                                     group_size_column = "group_size",
                                     min_group_size = 3) {
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c(group_size_column) %in% colnames(x)))
  group_size <- x[[ group_size_column ]]
  flt <- group_size >= min_group_size
  x[flt, ]
}


#' Remove features without good MS/MS spectra
#'
#' @param x MS-DIAL alignment results table
#' @param feat_group_column name of feature group column, default: "feat_group"
#' @param ms2_column name of MS/MS spectrum column, default: "ms2"
#' @param min_tic minimum summed intensity required for MS/MS spectra to be kept
#' @param min_ions minimum number of ions required for MS/MS spectra to be kept
#' @param keep_all_adducts if TRUE, keep all adducts of a feature group as long
#'   as one MS/MS spectrum fulfills above criteria. If FALSE (the default), keep
#'   only the adduct with the most intense MS/MS spectrum.
#'
#' @returns an object of the same class as 'x'
#' @export
filterAdductGroupsByMS2Quality <- function(x,
                                           feat_group_column = "feat_group",
                                           ms2_column = "ms2",
                                           min_tic = 1e3,
                                           min_ions = 3,
                                           keep_all_adducts = FALSE) {
  rank1 <- function(x) {
    rank(x, ties.method = "first")
  }
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c(feat_group_column, ms2_column) %in% colnames(x)))
  feat_group <- x[[ feat_group_column ]]
  ms2 <- x[[ ms2_column ]]
  s_ms2 <- prepareSpectra(x[, ms2_column], ms2_column)[[ "s" ]]
  tic_ms2 <- sapply(s_ms2, function(s) max(0, sum(s[,2]), na.rm = TRUE))
  nion_ms2 <- sapply(s_ms2, function(s) sum(is.finite(s[,2]) & s[,2] > 0))
  ms2_ok <- tic_ms2 >= min_tic & nion_ms2 >= min_ions
  grp_ok <- stats::ave(ms2_ok, feat_group, FUN = any)
  x_out <- if(keep_all_adducts) {
    x[grp_ok, ]
  } else {
    best_ms2 <- as.logical(stats::ave(-tic_ms2[grp_ok & ms2_ok], 
                                      feat_group[grp_ok & ms2_ok], 
                                      FUN = function(x) rank1(x) == 1))
    x[grp_ok & ms2_ok, ][best_ms2, ]
  }
  x_out
}


#' Subtract blank levels
#'
#' @param x MS-DIAL alignment results table
#' @param blank_idx index of blank samples; \code{NULL} to use \code{attr(x,
#'   "msdial_sam")$file_type == "Blank"}
#' @param fun function used to calculate blank levels
#'
#' @returns an object of the same class as 'x'
#' @export
subtractBlankLevels <- function(x, 
                                blank_idx = NULL, 
                                fun = c("median", "average")[1]) {
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
#' @param min_int minimum intensity that a feature must have
#' @param min_frac minimum fraction of samples (within each group if
#'   sample_group is defined accordingly) that contain the feature at an intensity >= min_int
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
    sum(x > min_int) / length(x) > min_frac
  )
  if(!is.null(dim(flt))) 
    flt <- apply(flt, 2, any)
  x[flt, ]
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
