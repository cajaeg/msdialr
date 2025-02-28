#' Extract chromatogram
#'
#' Wrapper for \code{\link[HiResTEC]{getMultipleBPC}()} accepting lists of
#' 'xcmsRaw' objects.
#' @param xraw 'xcmsRaw' object or a list of them
#' @param mz mass vector or NULL to extract the full mass range
#' @param mz_dev allowed m/z deviations, see
#'   \code{\link[HiResTEC]{getMultipleBPC}()} for details
#' @param rt retention time point in seconds or NULL to extract the full RT
#'   range
#' @param rt_dev RT window
#' @param zeroVal value to use for intensities <= 0 (typically NA or 0)
#' @param smooth window size for moving average smoother, 0 = no smoothing
#' @param EIC if 'TRUE' return a (sum-based) EIC instead of a (maximum-based) BPC
#'
#' @return A matrix with scan wise (rows) intensities for all requested m/z's
#'   (columns), or a list of such matrices.
#' @export
#'
#' @examples
#' mzml_files <- list.files(system.file("extdata", package = "msdialr"),
#' "GCQTOF_PAH.*\\.mzML", full.names = TRUE)
#' xraw <- loadRaws(mzml_files)
#'
#' # total ion chromatogram (TIC)
#' tic <- getChrom(xraw, EIC = TRUE) 
#' plotChrom(tic)
#'
#' # full-range base peak chromatogram (BPC) 
#' bpc <- getChrom(xraw, EIC = FALSE) 
#' plotChrom(bpc)
#'
#' # extracted ion chromatogram (XIC)
#' xic <- getChrom(xraw, mz = c(202.078, 101.042), mz_dev = 0.005, rt = 20*60, rt_dev = 1*60)
#' plotChrom(xic)
#'
#' par(mfrow = c(1,2))
#' plotChrom(xic, set.mfrow = FALSE) # override default multi-figure layout
#' 
#' # change scope for scaling
#' xic_scaled <- scaleChrom(xic, scope = "by_file")
#' plotChrom(xic_scaled)
getChrom <- function(xraw,
                   mz = NULL,
                   mz_dev = 0.01,
                   rt = NULL,
                   rt_dev = 2,
                   zeroVal = NA,
                   smooth = 0,
                   EIC = FALSE) {
  checkFinite <- function(x)
    if (!is.null(x) && any(!is.finite(x)))
      NULL
  else
    x
  mz <- checkFinite(mz)
  rt <- checkFinite(rt)
  if (inherits(xraw, "list")) {
    .n <- length(xraw)
    flt <- sapply(xraw, inherits, "xcmsRaw")
    out <- vector("list", .n)
    for (i in which(flt)) {
      out[[i]] <- getChrom(
        xraw[[i]],
        mz = mz,
        mz_dev = mz_dev,
        rt = rt,
        rt_dev = rt_dev,
        zeroVal = zeroVal,
        smooth = smooth,
        EIC = EIC
      )
    }
    return(out)
  } else {
    out <- HiResTEC::getMultipleBPC(
      xraw,
      mz = mz,
      mz_dev = mz_dev,
      rt = rt,
      rt_dev = rt_dev,
      zeroVal = zeroVal,
      smooth = smooth,
      returnEIC = EIC
    )
    if (is.null(out)) {
      out <- matrix(nrow = 0, ncol = length(mz))
      attr(out, "mz") <- mz
      attr(out, "mz_dev") <- mz_dev
    }
    attr(out, "from_file") <- basename(xraw@filepath[])
    return(out)
  }
}

#' @rdname getChrom
#' @export
getTIC <- function(xraw,
                   zeroVal = NA,
                   smooth = 0) {
  getChrom(xraw,
           zeroVal = zeroVal,
           smooth = smooth,
           EIC = FALSE)
}

#' @rdname getChrom
#' @param ... arguments passed to further functions
#' @export
plotTIC <- function(xraw,
                    zeroVal = NA,
                    smooth = 0,
                    ...) {
  tic <- getTIC(xraw, zeroVal = zeroVal, smooth = smooth)
  plotChrom(tic, ...)
}


#' Plot extracted ion chromatograms (EIC, BPC, TIC)
#'
#' @param chrom chromatogram as obtained by \code{\link{getChrom}()}
#' @param label_peaks label major peaks with retention times
#' @param fwhm full-width half maximum for peak detection (passed to
#'   \code{\link[xcms]{peaksWithMatchedFilter}()})
#' @param label.k number of chromatogram section to analyze: in each section the
#'   major peak will be labelled
#' @param legend display a legend of m/z values
#' @param max.legend how many m/z's to show in legend
#' @param set.mfrow set multi-figure layout according to number of chromatograms
#' @param expand.ylim include some extra space for peak labels
#' @param ... passed to underlying function \code{graphics::matplot()}
#'
#' @return (invisibly) original 'chrom' object with attribute 'peaks' attached if
#'   peak detection was performed
#' @export
#'
#' @examples
#' # see getChrom()
plotChrom <-
  function(chrom,
           label_peaks = TRUE,
           fwhm = 5,
           label.k = 5,
           legend = NULL,
           max.legend = 15,
           set.mfrow = TRUE,
           expand.ylim = 1.05,
           ...) {
    if (inherits(chrom, "list")) {
      if (set.mfrow) {
        opar <- graphics::par(mfrow = grDevices::n2mfrow(length(chrom)))
        on.exit(graphics::par(opar))
      }
      args <- list(...)
      argNames <- names(args)
      flt <- sapply(chrom, inherits, "matrix")
      gylim <- c(0, max(0, unlist(lapply(chrom[flt], as.vector)), na.rm = TRUE) * expand.ylim)
      args <- checkArgs(args, "ylim", gylim)
      args <- checkArgs(args, "xlab", "")
      args <- checkArgs(args, "ylab", "")
      args[["label_peaks"]] <- label_peaks
      args[["fwhm"]] <- fwhm
      args[["label.k"]] <- label.k
      args[["legend"]] <- legend
      args[["max.legend"]] <- max.legend
      args[["expand.ylim"]] <- expand.ylim
      for (i in 1:length(chrom)) {
        if (flt[i]) {
          args[["chrom"]] <- chrom[[i]]
          do.call(plotChrom, args)
          if (!is.null(attr(chrom[[i]], "from_file")))
            graphics::mtext(
              sub("(.+)(\\..+)$", "\\1", basename(attr(
                chrom[[i]], "from_file"
              ))),
              side = 3,
              adj = 0.01,
              cex = 0.66
            )
        } else {
          emptyplot("no data")
        }
      }
      return(invisible(NULL))
    }
    args <- list(...)
    argNames <- names(args)
    args <- checkArgs(args, "type", "l")
    args <- checkArgs(args, "lty", 1:5)
    args <- checkArgs(args, "col", 1:6)
    if (!is.null(.tmp <- attr(chrom, "rt"))) {
      chrom_rt <- .tmp / 60
      xlab <- "RT (min)"
    } else if ("x" %in% argNames) {
      chrom_rt <- args[["x"]] / 60
      xlab <- "RT (min)"
    } else {
      warning("no RT provided, using index")
      chrom_rt <- seq_len(nrow(chrom))
      fwhm <- 0.1
      xlab <- "Index"
    }
    args[["x"]] <- chrom_rt
    args <- checkArgs(args, "xlab", xlab)
    args[["y"]] <- chrom
    args <- checkArgs(args, "ylab", "Intensity (a.u.)")
    args <- checkArgs(args, "ylim", c(0, max(1, chrom, na.rm = TRUE) * expand.ylim))
    if (is.null(args[["ylim"]]))
      args[["ylim"]] <- c(0, max(1, chrom, na.rm = TRUE) * expand.ylim)
    if (nrow(chrom) > 0) {
      do.call(graphics::matplot, args = args)
      if (!any(is.finite(chrom)))
        graphics::text(stats::median(graphics::par("usr")[1:2]), 
                       stats::median(graphics::par("usr")[3:4]), 
                       labels = "no data")
      if (label_peaks) {
        rt <- chrom_rt * 60
        n_ions <- ncol(chrom)
        int <- if (n_ions == 1)
          as.vector(chrom[, 1])
        else
          apply(chrom, 1, function(x)
            max(c(0, x), na.rm = TRUE))
        pks <- data.frame(suppressWarnings(
          xcms::peaksWithMatchedFilter(
            int = int,
            rt = rt,
            fwhm = fwhm,
            max = 50,
            snthresh = 2
          )
        ))
        if (nrow(pks) > 0) {
          pks$rt <- pks$rt / 60
          pks <- pks[order(pks$rt), ]
          pks$grp <- hcgroup(pks$rt, k = min(nrow(pks), label.k))
          flt <- unlist(tapply(-pks$maxo, pks$grp, function(x)
            rank(x, ties.method = "first") == 1, simplify = FALSE))
          pks <- pks[flt, ]
          graphics::text(
            pks$rt,
            pks$maxo,
            labels = round(pks$rt, 2),
            pos = 3,
            offset = 0.1
          )
        }
      }
      show.legend <- (is.logical(legend) && legend) ||
        is.character(legend) || is.numeric(legend) ||
        (is.null(legend) && ncol(chrom) > 1)
      show.legend <- show.legend && !is.null(attr(chrom, "mz"))
      if (show.legend) {
        legtext <- sprintf("%.3f", round(attr(chrom, "mz"), 3))
        if (is.character(legend) || is.numeric(legend)) {
          legtext <- paste(legtext, legend)
        }
        if (!is.null(legtext)) {
          maxleg <- max.legend
          legtext <- legtext[1:min(maxleg, length(legtext))]
          graphics::legend(
            "topright",
            bty = "n",
            legend = legtext,
            lty = args[["lty"]],
            col = args[["col"]]
          )
        }
      }
    } else {
      emptyplot("no data")
    }
    rv <- chrom
    attr(rv, "rt") <- chrom_rt * 60
    attr(rv, "peaks") <- if (label_peaks &&
                             exists("pks", envir = environment()))
      pks
    else
      NULL
    invisible(rv)
  }



#' Scale chromatograms
#'
#' @param chrom chromatogram(s) as obtained by \code{\link{getChrom}()}
#' @param scope one of "global" (default), "by_file", "by_trace" or "none"
#'
#' @return matrix or list equivalent to original 'chrom'
#' @export
#'
#' @examples
#' # see getChrom()
scaleChrom <- function(chrom,
                       scope = c("global", "by_file", "by_trace", "none")[1]) {
  if (inherits(chrom, "list")) {
    scope <- match.arg(scope, choices = c("global", "by_file", "by_trace", "none"))
    return(lapply(chrom, scaleChrom, scope = switch(
      scope,
      global = max(1, unlist(chrom), na.rm = TRUE),
      by_file = "by_file",
      by_trace = "by_trace",
      none = "none"
    )))
  } else {
    old_attr <- attributes(chrom)
    out <- if (is.numeric(scope)) {
      chrom / scope * 100
    } else {
      scope <- match.arg(scope, choices = c("global", "by_file", "by_trace", "none"))
      switch(
        scope,
        global = chrom / max(1, chrom, na.rm = TRUE) * 100,
        by_file = chrom / max(1, chrom, na.rm = TRUE) * 100,
        by_trace = apply(chrom, 2, function(x)
          x / max(1, x, na.rm = TRUE) * 100),
        none = chrom
      )
    }
    mostattributes(out) <- old_attr
    return(out)
  }
}


#' Add extracted ion chromatograms to alignment results
#'
#' @param x data.frame as created by \code{\link{loadAlignmentResults}()}
#' @param in_column name of column containing m/z's to extract: either a list of
#'   spectra or single m/z values
#' @param rt_column name of column containing retention times
#' @param out_column name of results column
#' @param xraw list of 'xcmsRaw' objects. If NULL (the default), \code{attr(x,
#'   "raw_data")} is used. See \code{\link{addRaws}()}.
#' @param mz_dev width of m/z window to extract (passed to
#'   \code{\link{getChrom}()})
#' @param rt_dev width of RT window to extract (passed to
#'   \code{\link{getChrom}()})
#' @param zeroVal value to replace NA's by (passed to \code{\link{getChrom}()})
#' @param smooth (passed to \code{\link{getChrom}()})
#' @param max_mz restrict BPC extraction to the 'max_mz' highest peaks (if
#'   'in_column' contains spectra)
#' @param EIC if 'TRUE' add (sum-based) EICs instead of (maximum-based) BPCs
#' @param row_index index of rows to process, default NULL = process all rows. 
#' @param verbose display a progress bar, defaults to the output of
#'   \code{interactive()}
#' @param pivot_longer transform resulting BPC matrices to long format. Useful
#'   for plotting with 'ggplot'.
#' @return data.frame
#' @export
#'
#' @examples
#' # see ?loadAlignmentResults
#' 
addXICs <- function(x,
                    in_column = "s",
                    rt_column = c("average_rt_min", "rt_min", "rt"),
                    out_column = "xic",
                    xraw = NULL,
                    mz_dev = 0.01,
                    rt_dev = 10,
                    zeroVal = NA,
                    smooth = 0,
                    max_mz = 5,
                    EIC = FALSE,
                    row_index = NULL,
                    verbose = NULL,
                    pivot_longer = FALSE) {
  stopifnot(in_column %in% colnames(x))
  stopifnot(any(rt_column %in% colnames(x)))
  if (is.null(verbose))
    verbose <- interactive()
  if (is.null(xraw))
    xraw <- attr(x, "raw_data")
  stopifnot(!is.null(xraw))
  if (inherits(xraw, "xcmsRaw"))
    xraw <- list(xraw)
  if(is.null(row_index))
    row_index <- rep(TRUE, nrow(x))
  xmz <- x[[in_column]][row_index]
  mz <- if (is.list(xmz))
    # in_column is a list of spectra
    lapply(xmz, function(spec)
      spec[order(-spec[, 2]), ][1:min(max_mz, nrow(spec)), ][, 1])
  else # in_column is a numeric vector
    mz <- as.numeric(unlist(xmz))
  rt_column <- rt_column[which(rt_column %in% colnames(x))[1]]
  rt <- x[[rt_column]][row_index]
  out <- vector("list", length(mz))
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0,
                                max = length(out),
                                style = 3)
    on.exit(close(pb))
  }
  for (i in 1:length(out)) {
    if (verbose) {
      utils::setTxtProgressBar(pb, i)
    }
    tmp <- vector("list", length(xraw))
    flt <- sapply(xraw, inherits, "xcmsRaw")
    for (j in seq_along(xraw)[flt]) {
      tmp[[j]] <- getChrom(
        xraw[[j]],
        mz = mz[[i]],
        mz_dev = mz_dev,
        rt = rt[[i]] * 60,
        rt_dev = rt_dev,
        zeroVal = zeroVal,
        smooth = smooth,
        EIC = EIC
      )
    }
    names(tmp)[flt] <- sub("\\.mzML$", "", basename(sapply(xraw[flt], methods::slot, "filepath")))
    out[[i]] <- if (pivot_longer) {
      plyr::ldply(tmp, function(x)
        data.frame(rt = attr(x, "rt"), x, check.names = F), .id = "fromFile") |>
        tidyr::pivot_longer(
          cols = !tidyselect::matches("fromFile|rt"),
          names_to = "mz",
          values_to = "i"
        ) |>
        dplyr::mutate(mz = factor(round(as.numeric(mz), 4), levels = round(attr(tmp[[1]], "mz"), 4)))
    } else
      tmp
  }
  x[[out_column]] <- vector("list", nrow(x))
  x[[out_column]][row_index] <- out
  x |>
    relocateIntensityColumns()
}
