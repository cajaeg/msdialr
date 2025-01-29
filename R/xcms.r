#' Add MS raw data as attribute to alignment results data.frame
#'
#' @param x data.frame as created by \code{\link{loadAlignmentResults}()}
#' @param ... arguments passed to \code{\link{loadRaws}()} 
#'
#' @return 'x' with attribute 'raw_data' attached 
#' @export
#'
#' @examples
#' # see ?loadAlignmentResults()
addRaws <- function(x, 
                    ...) {
  out <- NULL
  attr_name <- "raw_data"
  msdial_sam <- attr(x, "msdial_sam")
  if(!is.null(msdial_sam)) {
    raw_files <- msdial_sam$raw_file
    if(any(!is.na(raw_files))) {
      out <- vector("list", nrow(msdial_sam))
      which.ok <- which(!is.na(raw_files))
      out[ which.ok ] <- loadRaws(raw_files[ which.ok ], ...)
    }
  }
  attr(x, attr_name) <- out
  return(x)
}


#' Load MS raw data
#'
#' Wrapper for \code{\link[xcms]{xcmsRaw}()} accepting multiple MS files and
#' optionally applying an intensity filter.
#' @param msfiles mz(X)ML files to load
#' @param profstep profile matrix step size parameter used by
#'   \code{\link[xcms]{xcmsRaw}()}
#' @param minabs absolute intensity threshold
#' @param ... further arguments passed to \code{\link[xcms]{xcmsRaw}()}
#'
#' @return list of \code{xcmsRaw} objects
#' @export
#'
#' @examples
#' mzml_files <- list.files(system.file("extdata", package = "msdialr"),
#' "GCQTOF_PAH.*\\.mzML", full.names = TRUE)
#' xraw <- loadRaws(mzml_files)
#' xraw
loadRaws <- function(msfiles,
                     profstep = 0,
                     minabs = 0,
                     ...) {
  plyr::llply(
    msfiles,
    loadXcmsRaw,
    profstep = profstep,
    minabs = minabs,
    ...,
    .parallel = FALSE
  )
}


#' Intensity filter for 'xcmsRaw'
#'
#' Internal function used by \code{\link{loadRaws}()}
#' @param xraw 'xcmsRaw' object
#' @param minabs intensity threshold
#' @param verbose print diagnostic message
#'
#' @return 'xcmsRaw'
#' @noRd
#'
#' @examples
#' # see \code{\link{loadRaws}()}
filterXcmsRaw <- function(xraw,
                          minabs = 0,
                          verbose = TRUE) {
  ## remove intensity values below threshold
  ## TODO: add time filters
  ## TODO: add relative threshold
  n_mz0 <- length(xraw@env$mz)
  n_scan0 <- length(xraw@scanindex)
  if (!is.null(minabs) &&
      minabs > 0 && any(xraw@env$intensity < minabs)) {
    scidx <- rep(1:length(xraw@scanindex), diff(c(xraw@scanindex, length(xraw@env$mz))))
    flt <- xraw@env$intensity >= minabs
    keep.scans <- unique(scidx[flt])
    scidx.new <- match(unique(scidx[flt]), scidx[flt]) - 1
    if (any(!flt)) {
      out <- xcms::deepCopy(xraw)
      out@env$mz <- xraw@env$mz[flt]
      out@env$intensity <- xraw@env$intensity[flt]
      out@scanindex <- as.integer(scidx.new)
      out@mzrange <- range(xraw@env$mz)
      out@scantime <- xraw@scantime[keep.scans]
      out@acquisitionNum <- xraw@acquisitionNum[keep.scans]
      out@scanrange <- as.integer(range(seq_along(out@scanindex)))
      tic <- tapply(out@env$intensity, rep(1:length(scidx.new), diff(c(
        scidx.new, length(out@env$mz)
      ))), FUN = sum)
      out@tic <- as.integer(round(tic))
      n_mz <- length(out@env$mz)
      n_scan <- length(out@scanindex)
      if (verbose)
        message(
          sprintf(
            "kept %d of %d m/z values (%.1f%%) and %d of %d scans (%.1f%%)",
            n_mz,
            n_mz0,
            round(n_mz / n_mz0 * 100, 1),
            n_scan,
            n_scan0,
            round(n_scan / n_scan0 * 100, 1)
          )
        )
      return(out)
    }
  } else {
    if (verbose)
      message("Leaving object unchanged")
    return(xraw)
  }
}

#' Load MS raw data (single file)
#'
#' Internal function used by \code{\link{loadRaws}()}
#' @param msfile mz(X)ML files to load
#' @param profstep profile matrix step size parameter used by \code{\link[xcms]{xcmsRaw}()}
#' @param minabs absolute intensity threshold
#' @param ... further arguments passed to \code{\link[xcms]{xcmsRaw}()}
#'
#' @return 'xcmsRaw'
#' @noRd
#'
#' @examples
#' # see \code{\link{loadRaws}()}
loadXcmsRaw <- function(msfile,
                        profstep = 0,
                        minabs = 0,
                        ...) {
  if (is.character(msfile) && file.exists(msfile)) {
    xraw <- suppressPackageStartupMessages(xcms::xcmsRaw(msfile, profstep = profstep, ...))
    return(if (is.null(minabs) ||
               minabs == 0)
      xraw
      else
        filterXcmsRaw(xraw, minabs = minabs))
  } else {
    return(NULL)
  }
}


#' Get scan data for a given retention time
#'
#' Wrapper for \code{\link[xcms]{getScan}()} accepting a retention time and a
#' list of \code{xcmsRaw}'s
#' @param xraw 'xcmsRaw' object or list of 'xcmsRaw' objects
#' @param rt retention time in seconds
#' @param ... passed to \code{\link[xcms]{getScan}()}
#'
#' @return data.frame
#' @export
#'
#' @examples
#'  mzml_files <- list.files(system.file("extdata", package = "msdialr"),
#' "GCQTOF_PAH.*\\.mzML", full.names = TRUE)
#' xraw <- loadRaws(mzml_files)
#' rt_min <- 21
#' scan <- getScanAtRT(xraw, rt_min * 60)
#' plot(scan[[1]], type = "h")
getScanAtRT <- function(xraw, rt, ...) {
  if (inherits(xraw, "list")) {
    lapply(xraw, getScanAtRT, rt = rt, ...)
  } else {
    scanidx <- which.min(abs(xraw@scantime - rt))
    out <- data.frame(xcms::getScan(xraw, scanidx, ...))
    colnames(out) <- c("mz", "i")
    out
  }
}
