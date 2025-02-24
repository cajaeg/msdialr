#' Plot mass spectrum
#'
#' @param x matrix, data.frame or list of n spectra (see Details)
#' @param label.cex peak label magnification factor, see \code{\link[graphics]{par}()}
#' @param label.col peak label colour
#' @param label.k number of peaks to label (k as in \code{\link[stats]{cutree}})
#' @param peak.col peak colour
#' @param losses matrix of neutral loss relationships as returned by \code{\link{getNeutralLosses}()}
#' @param losses.col line/text colour for neutral losses
#' @param mirror show mirror plot. This is the default for n = 2 spectra. Set to FALSE to force "normal" plot.
#' @param precursor.mz add indicator for precursor m/z / neutral mass
#' @param expand.ylim Add some space in upper plotting region for peak labels.
#' @param ... further arguments passed to \code{\link[graphics]{plot.default}()}
#'
#' @return NULL
#' @export
#'
#' @examples
#' spec <- str2spec("39:74 41:150 65:49 77:84 91:72 95:142 107:275 135:999 136:103 150:262")
#' plotSpec(spec)
plotSpec <- function(x,
                   label.cex = 1,
                   label.col = 1,
                   label.k = 5,
                   peak.col = 1,
                   losses = NULL,
                   losses.col = 8,
                   mirror = NULL,
                   precursor.mz = NULL,
                   expand.ylim = 1.05,
                   ...) {
  checkArgs <- function(args, name, val) {
    if (!name %in% names(args))
      args[[name]] <- val
    return(args)
  }
  createDefaultMzLabel <- function(pd, label.k) {
    maxPeaks = 500 # limit clustering to 'maxPeaks' highest peaks
    sel <- rank(-abs(pd$i), ties.method = "first") <= maxPeaks
    label.k <- min(label.k, with(pd[sel, ], lengths(split(id, id)))) # make sure k < n_peaks
    pd$mzgrp <- 0
    pd$mzgrp[sel] <- with(pd[sel, ], unsplit(lapply(
      split(mz, id), FUN = hcgroup, k = label.k
    ), id))
    pd$flt <- FALSE
    pd$flt[sel] <- with(pd[sel, ], unsplit(lapply(
      split(-abs(i), list(id, mzgrp)), rank, ties.method = "first"
    ), list(id, mzgrp))) == 1
    pd$mzgrp2 <- NA
    pd$mzgrp2[pd$flt] <- with(pd, unsplit(lapply(
      split(mz[flt], id[flt]), FUN = hcgroup, h = 5
    ), id[flt]))
    pd$flt2 <- FALSE
    pd$flt2[pd$flt] <- with(pd, unsplit(
      lapply(split(-abs(i)[flt], list(id[flt], mzgrp2[flt])), rank, ties.method = "first"),
      list(id[flt], mzgrp2[flt])
    )) == 1
    pd$label <- with(pd, ifelse(flt & flt2, round(mz, 3), label))
    return(pd)
  }
  spec_ok <- inherits(x, c("list", "character")) ||
    (inherits(x, c("data.frame", "matrix")) && nrow(x) > 0)
  if (!spec_ok) {
    emptyplot("no data")
    return(invisible(NULL))
  }
  x.in <- x
  if (inherits(x, c("data.frame", "matrix"))) {
    x <- list(x)
  }
  nspec <- length(x)
  ##
  ## process arguments
  ##
  args <- list(...)
  argNames <- names(args)
  args <- checkArgs(args, "xlab", "m/z")
  args <- checkArgs(args, "ylab", "Intensity")
  args <- checkArgs(args, "main", if (!is.null(.name <- attr(x[[1]], "name", exact = TRUE)))
    .name
    else
      "")
  ##
  ## main
  ##
  mirror <- (is.null(mirror) &&
               nspec == 2) || (!is.null(mirror) && mirror)
  args <- checkArgs(args, "xlim", range(unlist(lapply(x, function(x)
    x[, 1]))) * c(0.95, 1.05))
  gylim <- c(0, max(1, abs(unlist(
    lapply(x, function(x)
      x[, 2])
  )), na.rm = TRUE) * expand.ylim)
  if (nspec == 2 && mirror) {
    x[[2]][, 2] <- -x[[2]][, 2]
    gylim[1] <- -gylim[2]
  }
  args <- checkArgs(args, "ylim", gylim)
  if (nspec > 2 || (nspec == 2 && !mirror)) {
    opar <- graphics::par(mfrow = grDevices::n2mfrow(nspec))
    on.exit(graphics::par(opar))
    for (i in 1:length(x)) {
      args[["x"]] <- x[[i]]
      do.call(plotSpec, args = args)
      if (!is.null(attr(x[[i]], "from_file")))
        graphics::mtext(
          sub("(.+)(\\..+)$", "\\1", basename(attr(
            x[[i]], "from_file"
          ))),
          side = 3,
          adj = 0.01,
          ## line = -1,
          cex = 0.66
        )
    }
  } else {
    x <- lapply(x, function(obj) {
      if (!checkSpec(obj))
        asSpec(obj)
      else
        obj
    })
    pdat <- data.frame(do.call(rbind, x), id = rep(1:length(x), sapply(x, nrow)))
    pdat$peakcol[is.na(pdat$peakcol)] <- peak.col
    pdat$labelcol[is.na(pdat$labelcol)] <- label.col
    if (all(is.na(pdat$label)))
      pdat$label <- createDefaultMzLabel(pdat, label.k = min(label.k, nrow(pdat)))$label
    args[["x"]] <- pdat$mz
    args[["y"]] <- pdat$i
    args[["type"]] <- "h"
    args[["col"]] <- pdat$peakcol
    if (is.null(args[["ylim"]]))
      args[["ylim"]] <- c(0, max(1, pdat$i) * expand.ylim)
    args[["axes"]] <- FALSE
    do.call(graphics::plot.default, args = args)
    graphics::axis(1, at = graphics::axTicks(1), labels = graphics::axTicks(1))
    graphics::axis(2, at = graphics::axTicks(2), labels = abs(graphics::axTicks(2))) # don't use negative y values in mirror plot
    graphics::box()
    if (any(!is.na(pdat$label))) {
      with(
        pdat[!is.na(pdat$label), ],
        graphics::text(
          mz,
          i,
          labels = label,
          col = labelcol,
          pos = ifelse(sign(i) == 1, 3, 1),
          offset = 0.1,
          cex = label.cex
        )
      )
    }
    if (nspec == 1 &&
        !is.null(losses) &&
        ((is.matrix(losses) ||
          is.data.frame(losses)) && nrow(losses) >= 1)) {
      .x <- .y <- data.matrix(losses[, 1:2])
      .x[] <- pdat[, 1][.x]
      .y[] <- pdat[, 2][.y]
      graphics::segments(
        x0 = .x[, 1],
        x1 = .x[, 2],
        y0 = .y[, 1],
        y1 = .y[, 2],
        col = losses.col,
        lty = 3
      )
      graphics::text(
        apply(.x, 1, mean),
        apply(.y, 1, mean),
        labels = losses[, 3],
        col = losses.col,
        cex = 0.66
      )
    }
    if (!is.null(precursor.mz) && is.numeric(precursor.mz)) {
      graphics::points(
        x = precursor.mz,
        y = 0 - abs(diff(range(graphics::par(
          "usr"
        )[3:4]))) * 0.01,
        pch = 17,
        cex = 1.5
      )
    }
  }
}


newSpec <-
  function(mz = NULL,
           i = NULL,
           adduct = NULL,
           isogr = NULL,
           iso = NULL,
           charge = NULL,
           formula = NULL,
           exactmass = NULL,
           error = NULL,
           label = NULL,
           peakcol = NULL,
           labelcol = NULL,
           .nrow = NULL) {
    ## define and create a standard dataframe for mass spectra
    colDef <- c(
      "mz" = "double",
      "i" = "double",
      "adduct" = "character",
      "isogr" = "integer",
      "iso" = "integer",
      "charge" = "integer",
      "formula" = "character",
      "exactmass" = "double",
      "error" = "double",
      "label" = "character",
      "peakcol" = "character",
      "labelcol" = "character"
    )
    colNames <- names(colDef)
    nCol <- length(colDef)
    nRow <- if (is.null(.nrow))
      max(
        length(mz),
        length(i),
        length(isogr),
        length(iso),
        length(charge),
        length(formula),
        length(exactmass),
        length(error),
        length(label),
        length(peakcol),
        length(labelcol)
      )
    else
      as.integer(.nrow)
    out <- data.frame(matrix(nrow = nRow, ncol = nCol), stringsAsFactors = FALSE)
    colnames(out) <- colNames
    for (j in 1:length(colDef)) {
      x <- colNames[j]
      if (!is.null(get(x)))
        out[[x]] <- get(x)
      mode(out[[x]]) <- colDef[j]
    }
    return(out)
  }


#' Format mass spectrum
#'
#' @param x matrix or 
#'
#' @return data.frame
#' @export
#'
asSpec <- function(x) {
  stopifnot(inherits(x, c("matrix", "data.frame", "character")))
  MS <- newSpec() # get definition
  if (inherits(x, "character"))
    x <- str2spec(x)
  if (!inherits(x, "data.frame"))
    x <- data.frame(x)
  colnames(x)[1:2] <- colnames(MS)[1:2]
  newSpec(
    x$mz,
    x$i,
    x$adduct,
    x$isogr,
    x$iso,
    x$charge,
    x$formula,
    x$exactmass,
    x$error,
    x$label,
    x$peakcol,
    x$labelcol
  )
}


checkSpec <- function(x,
                      checkColnames = T,
                      checkMode = F) {
  colDef <- c(
    "mz" = "double",
    "i" = "double",
    "adduct" = "character",
    "isogr" = "integer",
    "iso" = "integer",
    "charge" = "integer",
    "formula" = "character",
    "exactmass" = "double",
    "error" = "double",
    "label" = "character",
    "peakcol" = "character",
    "labelcol" = "character"
  )
  ok <- inherits(x, "data.frame")
  if (ok && checkColnames) {
    ok <- all(names(colDef) %in% colnames(x))
  }
  if (ok && checkMode) {
    ok <- all(sapply(1:length(colDef), function(i)
      mode(x[[i]] == colDef[i])))
  }
  return(ok)
}


#' Detect neutral losses in mass spectrum
#'
#' @param x mass spectrum as matrix, data.frame or similar
#' @param nldef neutral loss definition. If NULL, uses 50 typical EI losses
#' @param mzabs allowed m/z deviation in amu
#' @param merdef mer definition for detecting multimers
#'
#' @return data.frame
#' @export 
#' 
getNeutralLosses <- function(x,
                             nldef = NULL,
                             mzabs = 0.005,
                             merdef = expand.grid(nmol = 2, dmz = 1.007276)) {
  if (is.null(nldef)) {
    neutrallosses <- NULL
    utils::data("neutrallosses", 
                envir = environment(), 
                package = "msdialr")
    nldef <- neutrallosses$EI
  }
  x.in <- x
  x <- data.matrix(x[, 1:2])
  mz <- x[, 1]
  nl_mass <- nldef$exactmass
  nl_name <- nldef$formula
  d <- as.matrix(stats::dist(mz))
  d[upper.tri(d, diag = T)] <- NA
  m <- abs(outer(d, nl_mass, "-")) < mzabs
  dimnames(m) <- list(round(mz, 4), round(mz, 4), nl_name)
  arr_idx <- plyr::ldply(apply(m, 3, which, arr.ind = T))
  if (nrow(arr_idx) > 0) {
    colnames(arr_idx) <- c("name", "mz1", "mz2")
  }
  else {
    arr_idx <- data.frame(matrix(
      nrow = 0,
      ncol = 3,
      dimnames = list(NULL, c("name", "mz1", "mz2"))
    ))
  }
  if (!is.null(merdef)) {
    mers <- arr_idx[0, ]
    for (i in 1:nrow(merdef)) {
      m <- outer((mz - merdef$dmz[i]), (mz - merdef$dmz[i]) *
                   merdef$nmol[i], FUN = "-")
      m[upper.tri(m, diag = T)] <- NA
      tmp <- which(abs(m) < mzabs, arr.ind = T)
      if (nrow(tmp) > 0) {
        tmp <- data.frame(tmp, name = sprintf("%dM_%d", merdef$nmol[i], round(merdef$dmz[i])))[, c(3, 1, 2)]
        colnames(tmp)[2:3] <- c("mz1", "mz2")
        mers <- rbind(mers, tmp)
      }
      else
        next
    }
    arr_idx <- rbind(arr_idx, mers)
  }
  if (nrow(arr_idx) > 0) {
    map <- which(abs(outer(x[, 1], x.in[, 1], "-")) <
                   .Machine$double.eps, arr.ind = TRUE)
    arr_idx[, 2:3] <- map[, 2][data.matrix(arr_idx[, 2:3])]
    res <- data.frame(arr_idx[, c(2, 3, 1)], type = "neutral_loss")
  }
  else {
    res <- data.frame(matrix(
      nrow = 0,
      ncol = 4,
      dimnames = list(NULL, c("mz1", "mz2", "name", "type"))
    ))
  }
  return(res)
}


#' Scale/normalize intensities in mass spectrum 
#'
#' @param x matrix or data.frame
#' @param newmax new maximum intensity (default: 100)
#' @param aggregationFun function to apply (default: max)
#'
#' @return matrix or data.frame
#' @export
#'
scaleSpec <- function(x,
                    newmax = 100,
                    aggregationFun = c("max", "sum")[1]) {
  aggregationFun <- match.arg(aggregationFun, c("max", "sum"))
  maxval <- eval(call(aggregationFun, x[, "i"], na.rm = T))
  x[, "i"] <- x[, "i"] / maxval * newmax
  return(x)
}

#' Clean mass spectrum
#'
#' @param x matrix or data.frame 
#' @param relthr relative intensity threshold
#' @param absthr absolute intensity threshold
#' @param dmz delta m/z
#' @param minmz minimum m/z 
#' @param maxmz maximum m/z
#'
#' @return matrix or data.frame
#' @export
#'
cleanSpec <- function(x,
                    relthr = 0,
                    absthr = 0,
                    dmz = NULL,
                    minmz = NULL,
                    maxmz = NULL) {
  x <- x[order(x[, "mz"]), , drop = FALSE]
  x <- x[!duplicated(x[, "mz"]), , drop = FALSE]
  x <- x[!is.na(x[, "mz"]), , drop = FALSE]
  x <- x[!is.na(x[, "i"]), , drop = FALSE]
  if (is.null(minmz))
    minmz <- min(x[, "mz"], na.rm = TRUE)
  if (is.null(maxmz))
    maxmz <- max(x[, "mz"], na.rm = TRUE)
  flt <-
    x[, "mz"] >= minmz &
    x[, "mz"] <= maxmz &
    x[, "i"] > (max(x[, "i"], na.rm = TRUE) * relthr) &
    x[, "i"] > absthr
  x <- x[flt, , drop = FALSE]
  if (!is.null(dmz) && nrow(x) > 1) {
    mzClus <- hcgroup(x[, "mz"], h = dmz)
    mzClusIntensityRank <- stats::ave(
      -x[, "i"],
      mzClus,
      FUN = function(x)
        rank(x, ties.method = "first")
    )
    flt <- mzClusIntensityRank == 1
    x <- x[flt, , drop = FALSE]
  }
  return(x)
}


#' Prepare spectra for annotation
#'
#' Convert spectra from character vector to data.frame representation, apply an intensity threshold and
#' add some columns to the data.frame for annotation.
#' @param x MS-DIAL results table as obtained by
#'   \code{\link{loadAlignmentResults}()} or \code{\link{loadConsoleResults}()}
#' @param in_column name of column to be converted, e.g. "ms1" (default), "ms2",
#'   or "spectrum"
#' @param out_column name of column containing the prepared spectra
#' @param intrel relative intensity threshold to be applied to each spectrum
#' @return data.frame
#' @export
#'
#' @examples
#' fp <- system.file("extdata/MSDconsole_single_result.msdial", package = "msdialr")
#' dat <- loadConsoleResults(fp)
#' colnames(dat) # spectrum column is 'spectrum'
#' dat <- prepareSpectra(dat, in_column = "spectrum")
#' head(dat$s[[1]])
#' plot(dat$s[[1]][,1:2], type = "h")
#'
#' # also see ?loadAlignmentResults
prepareSpectra <- function(x,
                           in_column = c("ms1", "spectrum"),
                           out_column = "s",
                           intrel = 0.001) {
  stopifnot(inherits(x, c("matrix", "data.frame", "tbl", "tbl_df")))
  stopifnot(any(in_column %in% colnames(x)))
  in_column <- in_column[which(in_column %in% colnames(x))[1]]
  s_chr <- x[[ in_column ]]
  s_mat <- lapply(s_chr, str2spec)
  s_mat <- lapply(s_mat, function(x) {
    x[, "i"] <- ifelse(is.finite(x[, "i"]), x[, "i"], 0)
    return(x)
  })
  s_mat <- lapply(s_mat, function(x) {
    x[x[, "i"] >= max(x[, "i"]) * intrel, ]
    return(x)
  })
  s_df <- lapply(s_mat, as.data.frame, stringsAsFactors = FALSE)
  s_df <- lapply(s_df, function(x) {
    x$isogr <- NA_integer_
    x$iso <- NA_integer_
    x$charge <- NA_integer_
    x$formula <- NA_character_
    x$exactmass <- NA_complex_
    x$error <- NA_character_
    x$label <- NA_character_
    x$peakcol <- 1
    x$labelcol <- 1
    return(x)
  })
  x[[ out_column ]] <- s_df
  x|>
    relocateIntensityColumns() |>
    updateIntensityColumnIndex()
}
