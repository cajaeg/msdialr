#' Import MS-DIAL alignment results
#'
#' Read "Height.txt" or "Area.txt" files as obtained by running "Export ->
#' Alignment result" in MS-DIAL. Supports file format created by MS-DIAL
#' versions >= 4.0.
#' @param x Alignment results file, e.g. "Height.txt", "Area.txt" or
#'   "AlignmentResult.msdial", or folder containing these files (in which case
#'   the last created file is used).
#' @param raw_folder Location of raw data. If given, file paths of mzML files
#'   will be attached to \code{attr(x, "msdial_sam")}. If NULL (the default),
#'   \code{dirname(txt)} will be checked for raw data. Pass "" to disable any
#'   search for raw data.
#' @param raw_suffix File type of raw data. Adjust this when file type differs
#'   from 'mzML'. Only relevant when 'raw_folder' is given. Search pattern is
#'   case-insensitive.
#' @param sample_names Optional character vector to use instead of original
#'   column names. No check for proper order is performed.
#' @param prefix Prefix added to alignment_id, e.g. "hilic_pos_". If NULL (the
#'   default), alignment_id's are kept unchanged.
#'
#' @return a data.frame with attribute "msdial_sam", a sample list derived from
#'   the file header containing sample class information as defined in MS-DIAL.
#'   This attribute can be manipulated e.g. to re-define sample groups or remove
#'   unwanted intensity columns. A copy of this attribute is attached as a
#'   second attribute "sam_orig", which is required to keep track of intensity
#'   and non-intensity columns. This attribute should not be manipulated.
#' 
#' @export
#' 
#' @examples
#' fp <- system.file("extdata", package = "msdialr")
#' height_file <- file.path(fp, "MSDIAL_Alignment_result_GCQTOF.txt")
#' msd <- loadAlignmentResults(height_file) |> dplyr::as_tibble()
#' head(msd)
#' head(getIntensityMatrix(msd)) # only the intensity matrix part
#' 
#' sam <- getSampleList(msd) # sample information as defined in MS-DIAL
#' 
#' sam1 <- sam[-2, ] # remove a sample
#' msd1 <- setSampleList(msd, sam1)
#' head(getIntensityMatrix(msd1)) # corresponding column is removed
#'
loadAlignmentResults <-
  function(x,
           raw_folder = NULL,
           raw_suffix = "mzml",
           sample_names = NULL,
           prefix = NULL
           ) {
    stopifnot(file.exists(x))
    if (file.info(x)$isdir) {
      fp <- x
      pat <- "Height.*[0-9]+\\.txt|Area.*[0-9]+\\.txt|AlignResult.*\\.msdial"
      fls <- findLatestFile(fp, pat)
      x <- fls
    }
    if (is.null(raw_folder))
      raw_folder <- dirname(x)
    hdr_nrow <- 4 # number of rows occupied by header
    ##
    ## evaluate header
    ##
    hdr <- utils::read.table(
      x,
      sep = "\t",
      header = F,
      skip = 0,
      as.is = T,
      check.names = F,
      comment.char = "",
      na.strings = c("NA", "No record", "null"),
      quote = "",
      nrows = hdr_nrow + 1
    )
    ok <- "Class" %in% hdr[1, ]
    if (!ok) {
      stop("\"Class\" keyword not found in header, check file format")
    }
    intensity_columns_start <- grep("^Class$", hdr[1, ])[1] + 1
    intensity_columns_end <- max(which(!is.na(hdr[2, ]) &
                                         !duplicated(unlist(hdr[5, ]))))
    intensity_columns <- intensity_columns_start:intensity_columns_end
    msdial_sam <- data.frame(t(hdr[, intensity_columns]), stringsAsFactors = FALSE)
    cn <- hdr[, grep("Class", hdr[1, ])[1]]
    cn <- gsub("[^a-zA-Z0-9]+", "_", tolower(trimws(cn)))
    cn <- sub("_$", "", cn)
    cn <- sub("^_", "", cn)
    cn[5] <- "id"
    colnames(msdial_sam) <- cn
    rownames(msdial_sam) <- NULL
    msdial_sam <- msdial_sam[, c(grep("^id$", colnames(msdial_sam)), which(!grepl("^id$", colnames(msdial_sam))))]
    other_columns <- 1:(intensity_columns[1] - 1)
    ##
    ## evaluate body
    ##
    out <- utils::read.table(
      x,
      sep = "\t",
      header = T,
      skip = hdr_nrow,
      as.is = T,
      check.names = F,
      comment.char = "",
      na.strings = c("NA", "No record", "null"),
      quote = ""
    )
    cn <- colnames(out)[other_columns]
    cn <- make.names(cn)
    cn <- gsub("[^a-zA-Z0-9]+", "_", tolower(trimws(cn))) # fix column names
    cn <- sub("_$", "", cn)
    cn <- sub("^_", "", cn)
    msd_mode <- if ("ei_spectrum" %in% cn)
      "GCMS"
    else
      "LCMS"
    if (msd_mode == "GCMS") {
      cn[cn == "ei_spectrum"] <- "ms1"
    } else {
      cn[cn == "ms1_isotopic_spectrum"] <- "ms1"
      cn[cn == "ms_ms_spectrum"] <- "ms2"
    }
    colnames(out)[other_columns] <- cn
    for (i in which(apply(out, 2, function(x)
      all(x %in% c("False", "True"))))) {
      out[, i] <- as.logical(out[, i]) # encode as logical where appropriate
    }
    if (!is.null(sample_names)) {
      if (length(sample_names) == length(intensity_columns)) {
        colnames(out)[intensity_columns] <- sample_names
        msdial_sam$id <- sample_names
      } else {
        stop("supplied sample names do not match number of sample columns")
      }
    }
    if (ncol(out) > max(intensity_columns)) {
      out <- out[, -seq(max(intensity_columns) + 1, ncol(out))] # remove Average & SD columns
    }
    if (raw_folder != "" && dir.exists(file.path(raw_folder))) {
      raw_suffix_pat <- sprintf("\\.%s$", raw_suffix)
      raw_files <- list.files(raw_folder,
                              raw_suffix_pat,
                              ignore.case = TRUE,
                              full.names = TRUE)
      pat <- sub(raw_suffix_pat, "", basename(raw_files), ignore.case = TRUE)
      o <- match(msdial_sam$id, pat)
      raw_files <- raw_files[o]
      msdial_sam$raw_file <- raw_files
    } else {
      msdial_sam$raw_file <- NA_character_
    }
    if (!is.null(prefix)) {
      out$alignment_id <- 
        formatNumericIDs(out$alignment_id, prefix = prefix)
    } 
    attr(out, "msdial_sam") <- msdial_sam
    attr(out, "sam_orig") <- msdial_sam
    message(sprintf(
      "%s data loaded for %d compounds x %d samples",
      msd_mode,
      nrow(out),
      length(intensity_columns)
    ))
    return(out)
  }


#' Get sample list
#' @rdname loadAlignmentResults
#'
#' @param x MS-DIAL alignment file
#'
#' @returns data.frame
#' @export
getSampleList <- function(x) {
  attr(x, "msdial_sam")
}


#' Set sample list
#' @rdname loadAlignmentResults
#'
#' @param x MS-DIAL alignment file
#' @param sam data.frame or similar, or NULL
#'
#' @returns 'x' with new or updated attribute 'msdial_sam'
#' @export
setSampleList <- function(x, sam) {
  stopifnot(inherits(sam, "data.frame"))
  attr(x, "msdial_sam") <- sam
  removeUnusedIntensityColumns(x, sam)
}


#' Get index of intensity columns
#' @rdname loadAlignmentResults
#' 
#' @param x MS-DIAL alignment file
#'
#' @return numeric()
#' @export
getIntensityColumns <- function(x) {
  id_column <- "id"
  sam <- getSampleList(x)
  out <- match(sam[[id_column]], colnames(x))
  if(any(!is.finite(out)))
    warning("sam not synchronized with colnames")
  out
}


#' Get intensity matrix part of MS-DIAL alignment results
#' @rdname loadAlignmentResults
#'
#' @param x MS-DIAL alignment file
#' @param as.matrix ensure that return value is a numeric matrix, default: TRUE
#'
#' @return a numeric matrix, or an object of the same class as 'x' for
#'   \code{as.matrix == FALSE}
#' @export
getIntensityMatrix <- function(x, as.matrix = TRUE) {
  cidx <- getIntensityColumns(x)
  if(as.matrix) {
    out <- data.matrix(x[, cidx])
    rownames(out) <- make.unique(as.character(x$alignment_id))
    out
  } else {
    x[, cidx]
  }
}


#' Move intensity columns to right edge of data frame
#' @rdname loadAlignmentResults
#'
#' @param x MS-DIAL alignment results table
#'
#' @returns an object of the same class as 'x'
#' @export
relocateIntensityColumns <- function(x) {
  col_idx_old <- getIntensityColumns(x)
  if(!is.null(col_idx_old)) {
    is_intensity_column <- (1:ncol(x)) %in% col_idx_old
    col_idx_new <- c(which(!is_intensity_column), which(is_intensity_column))
    x[, col_idx_new]
  } else {
    x
  }
}


#' Remove intensity columns not referenced in 'sam'
#'
#' @param x MS-DIAL alignment results table
#' @param sam sample table. If NULL, uses \code{attr(x, "msdial_sam")}
#'
#' @return an object of the same class as 'x'
#' @noRd
removeUnusedIntensityColumns <- function(x, sam = NULL) {
  id_column <- "id"
  if (is.null(sam)) {
    sam <- getSampleList(x)
  }
  if (!is.null(sam)) {
    sam_orig <- attr(x, "sam_orig")
    cidx <- match(sam_orig[[id_column]], colnames(x))
    cidx <- cidx[is.finite(cidx)]
    flt <- cidx[ !colnames(x)[cidx] %in% sam[[id_column]] ]
    if (length(flt) > 0) {
      message(sprintf("Removing %d of %d intensity columns", length(flt), length(cidx)))
      x[, -flt]
    } else {
      x
    }
  } else {
    x
  }
}


#' Create a new MSP library from MS-DIAL alignment results
#'
#' @param x alignment result table as created by loadAlignmentResult()
#' @param intabs intensity threshold applied to mass spectra (PEAKS field in the
#'   resulting msp)
#' @param make.unique apply a make.unique() to spectra names in case of
#'   duplicates, default TRUE
#' @param digits numeric(2) giving decimal digits to round m/z and intensity
#'   values to
#' @param ... (currently unused)
#' 
#' @return data.frame
#' @export
#'
#' @examples
#' fp <- system.file("extdata/MSDIAL_Alignment_result_GCQTOF.txt", package = "msdialr")
#' aligned <- loadAlignmentResults(fp)
#' flt <- aligned$metabolite_name != "Unknown" # restrict to identified compounds
#' aligned <- aligned[flt, ]
#' msp <- alignmentResult2msp(aligned)
#' \dontrun{
#' writeMSP(msp, file = "library.msp")
#' }
alignmentResult2msp <- function(x,
                                intabs = 100,
                                makeNamesUnique = TRUE,
                                digits = c(4, 0),
                                ...) {
  spectra <- lapply(x$ms1, function(x) {
    x <- str2spec(x)
    x <- if (any(x[, 2] >= intabs))
      x[x[, 2] >= intabs, , drop = FALSE]
    else
      x[which.max(x[, 2]), , drop = FALSE]
    x[, 1] <- round(x[, 1], digits = digits[1])
    x[, 2] <- round(x[, 2] / max(x[, 2], na.rm = TRUE) * 999, digits = digits[2])
    x
  })
  spectra_txt <- lapply(spectra, spec2str)
  out <- data.frame(matrix(nrow = nrow(x), ncol = 0))
  tmp <- x$metabolite_name
  if (any(duplicated(tmp)) && make.unique)
    tmp <- base::make.unique(tmp, sep = "_")
  out$Name <- tmp
  out$Formula <- x$formula
  out$RETENTIONTIME <- x$average_rt_min
  out$RETENTIONINDEX <- x$average_ri
  out$SMILES <- x$smiles
  out$INCHIKEY <- x$inchikey
  out$Ontology <- x$ontology
  out$`Num peaks` <- unlist(lapply(spectra, nrow))
  out$PEAKS <- unlist(spectra_txt)
  out
}


#' Remove text files created by an MS-DIAL run
#'
#' @param msdir MS-DIAL working directory
#' @param pat file pattern passed to \code{list.files()}. The default
#'   pattern selects .pai2, .dcl, .msdial and .aef files.
#' @param ask ask before deleting
#' @return (invisibly) character vector of file names
#' @export
#'
cleanMSDfiles <- function(msdir, pat = "\\.pai2$|\\.dcl$|\\.msdial$|\\.aef$", ask = TRUE) {
  fls <- list.files(msdir, pat, full.names = TRUE)
  if (length(fls) > 0) {
    print(fls)
    if (ask) {
      if (utils::askYesNo("Delete these files?")) {
        file.remove(fls)
      }
    } else {
      file.remove(fls)
    }
  } else {
    message("no files found")
  }
  invisible(fls)
}


#' Import alkane definition file
#' 
#' @param x text file with two tab-separated columns (Num|RT(min))
#' @param plot show diagnostic plot?
#' @return data.frame
#' @export
#'
#' @examples
#' fp <- system.file("extdata/Alkane.txt", package = "msdialr")
#' loadAlkanes(fp, plot = TRUE)
loadAlkanes <- function(x, plot = FALSE) {
  stopifnot(file.exists(x))
  if (file.info(x)$isdir) {
    fp <- x
    fls <- list.files(
      fp,
      pattern = "alkan.*\\.txt",
      full.names = TRUE,
      ignore.case = TRUE
    )
    stopifnot(length(fls) > 0)
    mtime <- file.info(fls)$mtime
    sel <- which(order(mtime, decreasing = TRUE) == 1)
    message(sprintf("Found %d alkane file(s) - using %s", 
                    length(fls), fls[sel]))
    x <- fls[sel]
  }
  out <- utils::read.csv(x, sep = "\t", as.is = TRUE)
  stopifnot(ncol(out) == 2)
  out$RI <- out[, 1] * 100
  colnames(out)[1:2] <- c("Num", "RT")
  if (plot)
    graphics::plot(out[, c("RI", "RT")])
  return(out)
}


#' Create an "MS-DIAL experiment information" file required for SWATH/MS^E data
#'
#' @param mzrange Method m/z range (scan range)
#' @param winwidth SWATH window width (Da)
#' @param ms2range m/z range to create SWATH windows for, often identical to
#'   \code{mzrange}
#' @param spectra_rate Only used to calculate accumulation time (\code{attr(out,
#'   "acc_time")})
#' @param ms1_acq_factor Increase accumulation time for MS1 scan by this factor
#' @param outfile If non-NULL, divert output to text file that can be used as
#'   MS-DIAL experiment information file
#'
#' @return data.frame, with additional attributes "cycle_time", "acc_time"
#'   (accumulation time) and "bruker_mrm_table". The latter can serve as
#'   template in Bruker DataAcquisition.
#' @export
#'
#' @examples
#' createSWATHtable()
#' \dontrun{
#' createSWATHtable(outfile = "MSDIAL_exp_info.txt")
#' }
#' out <- createSWATHtable(c(50, 600))
#' out
#' attr(out, "cycle_time")
#' attr(out, "acc_time")
#' attr(out, "bruker_mrm_table")
#' createSWATHtable(c(50, 500),
#'   winwidth = 25,
#'   ms1_acq_factor = 1.67,
#'   spectra_rate = 30) # reproduce HILIC SWATH scheme from Tsugawa et al. 2015 (cycle time 640 ms)
#' createSWATHtable(c(100, 1250),
#'   winwidth = 21,
#'   ms1_acq_factor = 10,
#'   spectra_rate = 90) # reproduce lipid SWATH scheme (cycle time 730 ms)
#' createSWATHtable(c(100, 1250),
#'   winwidth = 25,
#'   ms2range = c(350, 1250),
#'   ms1_acq_factor = 5,
#'   spectra_rate = 50) # compromise for Impact II (max scan rate=50)
createSWATHtable <- function(mzrange = c(50, 500),
                             winwidth = 25,
                             ms2range = NULL,
                             spectra_rate = 50,
                             ms1_acq_factor = 5,
                             outfile = NULL) {
  op <- options(stringsAsFactors = F)
  on.exit(options(op))
  if (is.null(ms2range))
    ms2range <- mzrange
  win_lower <- seq(from = ms2range[1],
                   to = ms2range[2] - winwidth,
                   by = winwidth)
  win_upper <- seq(from = ms2range[1] + winwidth,
                   to = ms2range[2],
                   by = winwidth)
  win_center <- seq(
    from = ms2range[1] + winwidth / 2,
    to = ms2range[2] - winwidth / 2,
    by = winwidth
  )
  nwin <- length(win_center)
  acq_factor <- c(ms1_acq_factor, rep(1, nwin))
  out_colnames <- c("Experiment", "MS Type", "Min m/z", "Max m/z")
  out <- data.frame(matrix(
    nrow = nwin + 1,
    ncol = length(out_colnames),
    dimnames = list(NULL, out_colnames)
  ),
  check.names = F)
  out[, 1] <- (1:nrow(out)) - 1
  out[, 2] <- c("SCAN", rep("SWATH", nwin))
  out[, 3] <- c(mzrange[1], win_lower)
  out[, 4] <- c(mzrange[2], win_upper)
  attr(out, "cycle_time") <- sum(acq_factor) / spectra_rate
  attr(out, "acc_time") <- acq_factor / spectra_rate
  attr(out, "bruker_mrm_table") <- data.frame(
    Mass = c(mzrange[1], win_center),
    Width = c(0, rep(winwidth, nwin)),
    Collision = c(8, rep(40, nwin)),
    xAcq = acq_factor
  )
  if (is.null(outfile)) {
    cat("Cycle time:", attr(out, "cycle_time"), "\n")
    cat("Accumulation times:", unique(attr(out, "acc_time")), "\n")
    return(out)
  } else {
    utils::write.table(
      out,
      file = outfile,
      sep = "\t",
      quote = F,
      col.names = T,
      row.names = F
    )
  }
}
