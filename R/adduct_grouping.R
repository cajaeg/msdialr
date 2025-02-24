#' Group co-eluting adduct peaks
#'
#' @param x MS-DIAL alignment results table
#' @param mz_column name of m/z column, default: "average_mz"
#' @param adduct_column name of adduct column, default: "adduct_type"
#' @param rt_column name of retention time column, default: "average_rt_min"
#' @param group_column column name to use for peak group, default: "adduct_group".
#'   Group labels follow scheme: <numeric ID>_<neutral mass>_<retention time>
#' @param delta_mz allowed m/z deviation (Da)
#' @param delta_rt allowed retention time deviation (min)
#' @param verbose display messages, defaults to the output of \code{interactive()}
#'
#' @returns an object identical to 'x', with additional column 'group_column'
#' @export
#' 
#' @examples
#' fp <- system.file("extdata/MSDIAL_Alignment_result_LC-MS.txt", package = "msdialr")
#' aligned <- loadAlignmentResults(fp)
#' length(unique(aligned$alignment_id)) # 230
#'
#' aligned1 <- aligned |>
#'   assignAdductGroups(verbose = TRUE)
#' length(unique(aligned1$adduct_group)) # 174
#' plotAdductGroups(aligned1)
#'
#' aligned2 <- aligned1 |>
#'   filterAdductGroupsBySize(min_group_size = 2, filter = TRUE)
#' length(unique(aligned2$nm_group)) # 36
#' plotAdductGroups(aligned2)
assignAdductGroups <- function(x, 
                               mz_column = "average_mz", 
                               adduct_column = "adduct_type",
                               rt_column = "average_rt_min",
                               group_column = "adduct_group",
                               delta_mz = 0.01,
                               delta_rt = 0.1, 
                               verbose = NULL) {
  if(is.null(verbose))
    verbose <- interactive()
  msg <- function(x, verbose) 
    if(verbose)
      message(x)
  length_unique <- function(x) 
    length(unique(x))
  mz_digits <- 3
  rt_digits <- 2
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c(mz_column, adduct_column, rt_column) %in% colnames(x)))
  x_mz <- x[[ mz_column ]]
  x_adduct_type <- x[[ adduct_column ]]
  x_rt <- x[[ rt_column ]]
  msg("Calculating neutral masses ...", verbose)
  nm <- mz2nm(x_mz, x_adduct_type)
  msg("Grouping neutral masses ...", verbose)
  nm_group <- hcgroup(nm, h = delta_mz)
  msg("Grouping retention times ...", verbose)
  rt_group <- hcgroup(x_rt, h = delta_rt)
  peak_group <- as.integer(factor(paste(nm_group, rt_group)))
  tbl <- table(peak_group)
  msg(sprintf("Found %d adduct groups containing between %d and %d adducts (median: %d)",
              length_unique(peak_group), 
              as.integer(min(tbl)), 
              as.integer(max(tbl)), 
              as.integer(median(tbl))), 
      verbose)
  nm_group_lbl <- round(stats::ave(nm, peak_group, FUN = stats::median), mz_digits)
  rt_group_lbl <- round(stats::ave(x_rt, peak_group, FUN = stats::median), rt_digits)
  peak_group_lbl <- sprintf(paste0("%s_%.", mz_digits, "f_%.", rt_digits, "f"), 
                            formatNumericIDs(peak_group), 
                            nm_group_lbl, 
                            rt_group_lbl)
  x[[ group_column ]] <- peak_group_lbl
  x |> 
    relocateIntensityColumns() |> 
    updateIntensityColumnIndex()
}


#' Plot adduct groups
#' @rdname assignAdductGroups
#'
#' @param x MS-DIAL alignment results table
#' @param mz_column name of m/z column, default "average_mz"
#' @param rt_column name of retention time column, default "average_rt_min"
#' @param adduct_column name of adduct column, default "adduct_type"
#'
#' @returns NULL
#' @export
plotAdductGroups <- function(x, 
                             mz_column = "average_mz",
                             rt_column = "average_rt_min",
                             adduct_column = "adduct_type") {
  x_mz <- x[[ mz_column ]]
  x_rt <- x[[ rt_column ]]
  x_adduct_type <- x[[ adduct_column ]]
  fct <- factor(x_adduct_type)
  tmp <- sort(table(fct), decreasing = TRUE)
  fct <- factor(fct, levels = levels(fct)[match(names(tmp), levels(fct))])
  bg_col <- match(fct, levels(fct))
  plot(x_rt, x_mz, type = "n", 
       xlab = "Retention Time (min)", ylab = "m/z")
  points(x_rt, x_mz, 
         col = NULL, pch = 21, bg = bg_col)
  legend("topleft", 
         legend = sprintf("%s (%d)", levels(fct), as.numeric(tmp)), 
         pch = 21, col = NA, pt.bg = seq_along(levels(fct)))
}


#' Flag/remove peak groups based on number of distinct adducts
#' @rdname assignAdductGroups
#'
#' @param x MS-DIAL alignment results table
#' @param group_column name of peak group column, default: "nm_group"
#' @param adduct_column name of adduct column, default: "adduct_type"
#' @param min_group_size minimum number of distinct adducts
#' @param flag_column name of column to use for flag, default: "passes_groupsize_filter"
#' @param filter filter results according to flag, default: FALSE
#' @param verbose display messages, defaults to the output of \code{interactive()}
#'
#' @returns an object of the same class as 'x', with additional column 'flag_column'
#' @export
filterAdductGroupsBySize <- function(x,
                                     group_column = "adduct_group",
                                     adduct_column = "adduct_type",
                                     min_group_size = 3,
                                     flag_column = "passes_groupsize_filter",
                                     filter = FALSE,
                                     verbose = NULL) {
  if (is.null(verbose))
    verbose <- interactive()
  msg <- function(x, verbose)
    if (verbose)
      message(x)
  length_unique <- function(x)
    length(unique(x))
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c(group_column, adduct_column) %in% colnames(x)))
  x_adduct_group <- x[[group_column]]
  x_adduct_type <- x[[adduct_column]]
  group_size <- as.integer(stats::ave(x_adduct_type, x_adduct_group, FUN = length_unique))
  flag <- group_size >= min_group_size
  x[[flag_column]] <- flag
  flt <- flag
  x_out <- if (filter)
    x[flt, ]
  else
    x
  msg(
    sprintf(
      "Found %d of %d adduct groups passing group size filter",
      length_unique(x_adduct_group[flt]),
      length_unique(x_adduct_group)
    ),
    verbose
  )
  x_out |>
    relocateIntensityColumns() |>
    updateIntensityColumnIndex()
}


#' Rank/remove adduct peaks based on MS/MS intensity
#' @rdname assignAdductGroups
#'
#' @param x MS-DIAL alignment results table
#' @param group_column name of feature group column, default: "adduct_group"
#' @param ms2_column name of MS/MS spectrum column, default: "ms2"
#' @param min_tic minimum summed intensity required for MS/MS spectra to be kept
#' @param min_ions minimum number of ions required for MS/MS spectra to be kept
#' @param rank_column name of column to use for rank, default: "rank_MS2".
#'   Spectra not passing above criteria are assigned \code{NA}, so filtering for
#'   \code{is.finite(rank_MS2) & rank_MS2 == 1} yields a clean table for
#'   further processing.
#' @param filter filter results according to rank, default: FALSE
#' @param verbose display messages, defaults to the output of
#'   \code{interactive()}
#'
#' @returns an object of the same class as 'x', with additional column
#'   'rank_column'
#' @export
filterAdductGroupsByMS2 <- function(x,
                                    group_column = "adduct_group",
                                    ms2_column = "ms2",
                                    min_tic = 1e3,
                                    min_ions = 3,
                                    rank_column = "rank_MS2",
                                    filter = FALSE,
                                    verbose = NULL) {
  if (is.null(verbose))
    verbose <- interactive()
  msg <- function(x, verbose)
    if (verbose)
      message(x)
  rank1 <- function(x) {
    rank(x, ties.method = "first")
  }
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c(group_column, ms2_column) %in% colnames(x)))
  x_adduct_group <- x[[group_column]]
  x_ms2 <- x[, ms2_column] |> 
    setSampleList(NULL) # avoid warning in 'prepareSpectra()'
  s_ms2 <- prepareSpectra(x_ms2, ms2_column)[["s"]]
  tic_ms2 <- sapply(s_ms2, function(s)
    max(0, sum(s[, 2]), na.rm = TRUE))
  nion_ms2 <- sapply(s_ms2, function(s)
    sum(is.finite(s[, 2]) & s[, 2] > 0))
  ms2_ok <- tic_ms2 >= min_tic & nion_ms2 >= min_ions
  ms2_rnk <- rep(NA_integer_, nrow(x))
  ms2_rnk[ms2_ok] <- as.integer(stats::ave(-tic_ms2[ms2_ok], x_adduct_group[ms2_ok], FUN = rank1))
  x[[ rank_column ]] <- ms2_rnk
  x_out <- if (filter) {
    x[is.finite(ms2_rnk) & ms2_rnk == 1, ]
  } else {
    x
  }
  n_ok <- length(unique(x_adduct_group[is.finite(ms2_rnk) &
                                         ms2_rnk == 1]))
  msg(sprintf(
    "Found %d of %d adduct groups passing MS2 filter",
    n_ok,
    length(unique(x_adduct_group))
  ),
  verbose)
  x_out |>
    relocateIntensityColumns() |>
    updateIntensityColumnIndex()
}

