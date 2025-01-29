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
#'   \code{dirname(txt)} will be checked for raw data. Pass "" to prevent any
#'   search for raw data.
#' @param raw_suffix File type of raw data. Adjust this when file type differs
#'   from 'mzML'. Only relevant when 'raw_folder' is given.
#' @param sample_names Optional character vector to use instead of original
#'   column names. No check for proper order is performed.
#'
#' @return data.frame with two attributes attached: "msdial_sam" and
#'   "intensity_columns".
#'
#'   \code{attr(x, "msdial_sam")} is derived from the file header and contains
#'   sample class information as defined in MS-DIAL.
#'
#'   \code{attr(x, "intensity_columns")} is an index of columns holding MS
#'   intensities (sample columns).
#' @export
#'
#' @examples
#' fp <- system.file("extdata", package = "msdialr")
#' height_file <- file.path(fp, "MSDIAL_Alignment_result_GCQTOF.txt")
#' msd <- loadAlignmentResults(height_file)
#' head(msd)
#' 
#' int_mat <- getIntensityMatrix(msd)
#' hist(apply(int_mat, 1, median))
loadAlignmentResults <-
  function(x,
           raw_folder = NULL,
           raw_suffix = "mzml",
           sample_names = NULL) {
    stopifnot(file.exists(x))
    if (file.info(x)$isdir) {
      fp <- x
      fls <- list.files(
        fp,
        pattern = "Height.*[0-9]+\\.txt|Area.*[0-9]+\\.txt|AlignResult.*\\.msdial",
        full.names = TRUE,
        ignore.case = TRUE
      )
      stopifnot(length(fls) > 0)
      mtime <- file.info(fls)$mtime
      sel <- which(order(mtime, decreasing = TRUE) == 1)
      message(sprintf("Found %d alignment results file(s) - using %s", 
                      length(fls), fls[sel]))
      x <- fls[sel]
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
    cn[5] <- "name"
    colnames(msdial_sam) <- cn
    rownames(msdial_sam) <- NULL
    msdial_sam <- msdial_sam[, c(grep("^name$", colnames(msdial_sam)), which(!grepl("^name$", colnames(msdial_sam))))]
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
    ## msdial_sam$name <- colnames(out)[intensity_columns] # include sample column names in 'msdial_sam'
    ## msdial_sam <- msdial_sam[,c(5,1:4)]
    for (i in which(apply(out, 2, function(x)
      all(x %in% c("False", "True"))))) {
      out[, i] <- as.logical(out[, i]) # encode as logical where appropriate
    }
    if (!is.null(sample_names)) {
      if (length(sample_names) == length(intensity_columns)) {
        colnames(out)[intensity_columns] <- sample_names
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
      o <- match(msdial_sam$name, pat)
      raw_files <- raw_files[o]
      msdial_sam$raw_file <- raw_files
    } else {
      msdial_sam$raw_file <- NA_character_
    }
    attr(out, "msdial_sam") <- msdial_sam
    attr(out, "intensity_columns") <- intensity_columns
    message(sprintf(
      "%s data loaded for %d compounds x %d samples",
      msd_mode,
      nrow(out),
      length(intensity_columns)
    ))
    return(out)
  }


#' Get index of intensity columns
#'
#' @param x MS-DIAL alignment file
#'
#' @return numeric()
#' @export
getIntensityColumns <- function(x) {
  attr(x, "intensity_columns")
}


#' Get intensity matrix portion of MS-DIAL alignment file
#'
#' @param x MS-DIAL alignment file
#'
#' @return numeric()
#' @export
getIntensityMatrix <- function(x) {
  x[, getIntensityColumns(x)]
}


# #' Update sample column index
# #'
# #' @param x MS-DIAL alignment file
# #' @param ref reference vector of sample names, typically from a separate sample list
# #'
# #' @return 'x' with updated attribute 'intensity_columns'
# #' @export
# updateSampleColumnIndex <- function(x, ref = NULL) {
#   if (!is.null(ref)) {
#     map <- match(ref, colnames(x))
#     if (any(!is.finite(map))) {
#       message(sprintf(
#         "%d sample(s) from 'ref' not found in 'x': %s\nthis/these should be removed from 'ref' to preserve index",
#         sum(!is.finite(map)),
#         paste(ref[!is.finite(map)], sep = "", collapse = ", ")
#       ))
#     }
#     icol_old <- getIntensityColumns(x)
#     if (length(icol_old) > sum(is.finite(map))) {
#       message(sprintf(
#         "%d sample(s) from 'x' not found in 'ref': %s\nremoving this/these from index",
#         sum(!colnames(x)[icol_old] %in% ref),
#         paste(colnames(x)[icol_old][!colnames(x)[icol_old] %in% ref],
#               sep = "", collapse = ", ")
#       ))
#     }
#     attr(x, "intensity_columns") <- map[is.finite(map)]
#     x
#   } else {
#     x
#   }
# }
# 


#' Create a new MSP library from MS-DIAL alignment results
#'
#' @param x alignment result table as created by loadAlignmentResult()
#' @param intabs intensity threshold applied to mass spectra (PEAKS field in the
#'   resulting msp)
#' @param makeNamesUnique apply a make.unique() to the spectra names
#'   (default)
#' @param digits numeric(2) giving decimal digits to round m/z and intensity values to
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
  if (any(duplicated(tmp)) && makeNamesUnique)
    tmp <- make.unique(tmp, sep = "_")
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


#' Import MS-DIAL console results
#'
#' Import '.msdial' files created by MS-DIAL console app.
#' @param x .msdial file for individual MS files. Use
#'   \code{\link{loadAlignmentResults}()} to import alignment results for
#'   multiple files.
#' @return data.frame
#' @export
#'
#' @examples
#' fp <- system.file("extdata/MSDconsole_single_result.msdial", package = "msdialr")
#' dat <- loadConsoleResults(fp)
loadConsoleResults <- function(x) {
  out <- utils::read.table(
    x,
    header = TRUE,
    quote = "",
    sep = "\t",
    check.names = FALSE,
    as.is = TRUE,
    comment.char = ""
  )
  cn <- colnames(out)
  cn <- gsub("[^a-zA-Z0-9]+", "_", trimws(tolower(cn)))
  cn <- sub("_+$", "", cn)
  cn <- sub("^_+", "", cn)
  colnames(out) <- cn
  peakid <- if (nrow(out) > 0)
    1:nrow(out)
  else
    numeric()
  out <- data.frame(peakid = peakid, out, stringsAsFactors = FALSE)
  message(sprintf("peak data loaded for %d compounds", nrow(out)))
  return(out)
}


#' Create a method file required by MS-DIAL console
#'
#' MS-DIAL parameters in console are set via a text file (option -m). Use this function to create that text file.
#' @param outfile name of text file to be created or 'NULL' to return method contents as list
#' @param data_type "Centroid" or "Profile"
#' @param ion_mode "Positive" or "Negative"
#' @param accuracy_type "IsAccurate" or "IsNominal"
#' @param retention_time_begin RT [min]
#' @param retention_time_end RT [min]
#' @param mass_range_begin m/z
#' @param mass_range_end m/z
#' @param number_of_threads 1
#' @param smoothing_method "LinearWeightedMovingAverage"
#' @param smoothing_level 2
#' @param average_peak_width 30
#' @param minimum_peak_height 1000
#' @param mass_slice_width 0.1 (for accurate mass)
#' @param mass_accuracy 0.025 (for accurate mass)
#' @param sigma_window_value 0.5
#' @param amplitude_cut_off 10
#' @param msp_file MSP file
#' @param ri_index_file ""
#' @param retention_type "RT" or "RI"
#' @param ri_compound "Alkanes" or "Fames"
#' @param retention_time_tolerance_for_identification 0.5
#' @param retention_index_tolerance_for_identification 20
#' @param ei_similarity_tolerance_for_identification 70
#' @param identification_score_cut_off 70
#' @param use_retention_information_for_scoring "True"
#' @param use_retention_information_for_filtering "False"
#' @param reference_file reference file for alignment
#' @param retention_time "RT" or "RI"
#' @param retention_time_tolerance_for_alignment 0.075
#' @param retention_index_tolerance_for_alignment 20
#' @param ei_similarity_tolerance 70
#' @param retention_time_factor 0.5
#' @param ei_similarity_factor 0.5
#' @param peak_count_filter 0
#' @param qc_at_least_filter "True"
#'
#' @return list or NULL depending on 'outfile'
#' @export
#'
#' @examples
#' createMSDCmethod() # print method contents to stdout
#' \dontrun{
#' createMSDCmethod(outfile = "test.txt") # write tag=value pairs
#' }
createMSDCmethod <- function(outfile = NULL,
                             ## Data type
                             data_type = "Centroid",
                             ion_mode = c("Positive", "Negative")[1],
                             accuracy_type = c("IsAccurate", "IsNominal")[1],
                             ## Data collection
                             retention_time_begin = 0,
                             retention_time_end = 100,
                             mass_range_begin = 0,
                             mass_range_end = 1000,
                             ## Data processing
                             number_of_threads = 1,
                             ## Peak detection
                             smoothing_method = "LinearWeightedMovingAverage",
                             smoothing_level = 2,
                             average_peak_width = 30,
                             minimum_peak_height = 1000,
                             mass_slice_width = 0.1,
                             mass_accuracy = 0.025,
                             ## MS1Dec parameters
                             sigma_window_value = 0.5,
                             amplitude_cut_off = 100,
                             ## Identification
                             msp_file = "",
                             ri_index_file = "",
                             retention_type = c("RT", "RI")[1],
                             ri_compound = c("Alkanes", "Fames")[1],
                             retention_time_tolerance_for_identification = 0.5,
                             retention_index_tolerance_for_identification = 20,
                             ei_similarity_tolerance_for_identification = 70,
                             identification_score_cut_off = 70,
                             use_retention_information_for_scoring = "True",
                             use_retention_information_for_filtering = "False",
                             ## Alignment
                             reference_file = "",
                             retention_time = c("RT", "RI")[1],
                             retention_time_tolerance_for_alignment = 0.15,
                             retention_index_tolerance_for_alignment = 200,
                             ei_similarity_tolerance = 70,
                             retention_time_factor = 0.5,
                             ei_similarity_factor = 0.5,
                             ## Filtering
                             peak_count_filter = 0,
                             qc_at_least_filter = "True") {
  method_contents <- list(
    ## Data type
    `Data type` = data_type,
    `Ion mode` = ion_mode,
    `Accuracy type` = accuracy_type,
    ## Data collection
    `Retention time begin` = retention_time_begin,
    `Retention time end` = retention_time_end,
    `Mass range begin` = mass_range_begin,
    `Mass range end` = mass_range_end,
    ## Data processing
    `Number of threads` = number_of_threads,
    ## Peak detection
    `Smoothing method` = smoothing_method,
    `Smoothing level` = smoothing_level,
    `Average peak width` = average_peak_width,
    `Minimum peak height` = minimum_peak_height,
    `Mass slice width` = mass_slice_width,
    `Mass accuracy` = mass_accuracy,
    ## MS1Dec parameters
    `Sigma window value` = sigma_window_value,
    `Amplitude cut off` = amplitude_cut_off,
    ## Identification
    `MSP file` = msp_file,
    `RI index file pathes` = if (identical(ri_index_file, ""))
      ""
    else
      "RI_index_file_paths.txt",
    `Retention type` = retention_type,
    `RI compound` = ri_compound,
    `Retention index tolerance for identification` = retention_index_tolerance_for_identification,
    `Retention time tolerance for identification` = retention_time_tolerance_for_identification,
    `EI similarity tolerance for identification` = ei_similarity_tolerance_for_identification,
    `Identification score cut off` = identification_score_cut_off,
    `Use retention information for scoring` = use_retention_information_for_scoring,
    `Use retention information for filtering` = use_retention_information_for_filtering,
    ## Alignment
    `Reference file` = reference_file,
    `Retention time` = retention_time,
    `Retention time tolerance for alignment` = retention_time_tolerance_for_alignment,
    `Retention index tolerance for alignment` = retention_index_tolerance_for_alignment,
    `EI similarity tolerance` = ei_similarity_tolerance,
    `Retention time factor` = retention_time_factor,
    `EI similarity factor` = ei_similarity_factor,
    ## Filtering
    `Peak count filter` = peak_count_filter,
    `QC at least filter` = qc_at_least_filter
  )
  if (!identical(ri_index_file, "")) {
    attr(method_contents, "ri_index_file_paths") <- "RI_index_file_paths.txt"
    attr(method_contents, "ri_index_file") <- ri_index_file
  }
  if (is.null(outfile)) {
    return(method_contents)
  } else {
    cat(sprintf("%s: %s", names(method_contents), unlist(method_contents)),
        sep = "\n",
        file = outfile)
    invisible(NULL)
  }
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


#' Trigger MS-DIAL console process
#'
#' Add \code{options(msdialr_msdialconsole_path = "/path/to/MSDIAL
#' ver.4.9.221218 Windowsx64/MsdialConsoleApp.exe")} (adjusted accordingly to your system) to your .Rprofile file. Use
#' \code{Sys.getenv("R_USER")} to check where this file is located on your
#' system.
#' @param msdir directory holding MS files supported by MS-DIAL (.abf, .mzML,
#'   .mzXML, .cdf etc.)
#' @param analysis_type one of "gcms" (default), "lcmsdda", "lcmsdia",
#'   "lcimmsdda", "lcimmsdia"
#' @param method list of parameters as returned by
#'   \code{\link{createMSDCmethod}()} or, alternatively, name of existing
#'   MS-DIAL method file (.txt).
#' @param readResults import results files after processing
#' @param skipRun don't run MS-DIAL console, just read existing .msdial files
#' @param verbose show output of MS-DIAL console
#' @return List of data.frame's containing MS-DIAL results, or a vector of file
#'   names for 'readResults' == FALSE
#' @export
#'
#' @examples
#' \dontrun{
#' mzml_files <- list.files(system.file("extdata", package = "msdialr"),
#'   "GCQTOF_PAH.*\\.mzML",
#'   full.names = TRUE)
#' workdir <- tempfile()
#' dir.create(workdir)
#' file.copy(mzml_files, workdir)
#' full_results <- runMSDconsole(workdir)
#' alignment_results <- full_results[[1]]
#' # check ?loadAlignmentsResults for examples of further processing these results
#' cleanMSDfiles(workdir, ask = FALSE)
#' }
#'
#' ## use a custom method file
#' \dontrun{
#' msdc_method <- createMSDCmethod(minimum_peak_height = 10000)
#' (result_files <- runMSDconsole(workdir, method = msdc_method, readResults = FALSE)
#' peaks <- loadConsoleResults(resfiles[[2]])
#' }
#'
#' ## clean up temporary folder
#' \dontrun{
#' unlink(workdir, recursive = TRUE, force = TRUE)
#' }
runMSDconsole <- function(msdir,
                          analysis_type = c("gcms", "lcmsdda", "lcmsdia", "lcimmsdda", "lcimmsdia")[1],
                          method = createMSDCmethod(),
                          readResults = TRUE,
                          skipRun = FALSE,
                          verbose = TRUE) {
  MSDCexe <- normalizePath(getOption("msdialr_msdialconsole_path"))
  if (is.null(MSDCexe) || !file.exists(MSDCexe))
    stop(
      "Set option 'msdialr_msdialconsole_path' to the full path of MSD console executable."
    )
  msdir <- normalizePath(msdir)
  stopifnot(dir.exists(msdir))
  if (!skipRun) {
    if (is.list(method)) {
      if (!is.null(attr(method, "ri_index_file_paths"))) {
        mzml_files <- normalizePath(list.files(
          msdir,
          "\\.mzml",
          full.names = TRUE,
          ignore.case = TRUE
        ))
        ri_index_file_path_outfile <- normalizePath(file.path(msdir, basename(
          attr(method, "ri_index_file_paths")
        )))
        cat(sprintf("%s\t%s", mzml_files, attr(method, "ri_index_file")), file = ri_index_file_path_outfile)
        method[["RI index file pathes"]] <- ri_index_file_path_outfile
      }
      method_file <- file.path(msdir, "MSDCparams.txt")
      cat(sprintf("%s: %s", names(method), unlist(method)),
          sep = "\n",
          file = method_file)
      
    } else if (is.character(method) && file.exists(method)) {
      method_file <- method
    }
    method_file <- normalizePath(method_file)
    stopifnot(file.exists(method_file))
    system2(
      MSDCexe,
      args = c(
        analysis_type,
        "-i",
        utils::shortPathName(msdir),
        "-o",
        utils::shortPathName(msdir),
        "-m",
        utils::shortPathName(method_file)
      ),
      stdout = if (verbose)
        ""
      else
        FALSE,
      stderr = if (verbose)
        ""
      else
        FALSE
    )
  }
  msdial_files <- list.files(msdir, "\\.msdial$", full.names = TRUE)
  if (readResults) {
    tryCatch({
      out <- vector("list", 2)
      names(out) <- c("aligned", "peaks")
      alignment_file <- grepl("AlignResult", msdial_files)
      if (any(alignment_file)) {
        fls <- msdial_files[which(alignment_file)]
        alignment_results <- lapply(fls, loadAlignmentResults, raw_folder = msdir)
        names(alignment_results) <- basename(sub("\\.msdial$", "", fls))
        out[[1]] <- alignment_results
      }
      single_files <- !alignment_file
      if (any(single_files)) {
        fls <- msdial_files[which(single_files)]
        single_results <- lapply(fls, loadConsoleResults)
        names(single_results) <- basename(sub("\\.msdial$", "", fls))
        out[[2]] <- single_results
      }
      return(out)
    }, error =
      function(e) {
        message("Could not read results, return file names only")
        return(msdial_files)
      })
  } else {
    return(msdial_files)
  }
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


#' Move intensity columns to right edge of data frame
#'
#' @param x MS-DIAL alignment results table
#'
#' @returns an object of the same class as 'x'
#' @noRd
relocateIntensityColumns <- function(x) {
  attr_intensity_columns <- "intensity_columns"
  col_idx_old <- attr(x, attr_intensity_columns)
  is_intensity_column <- (1:ncol(x)) %in% col_idx_old
  col_idx_new <- c(which(!is_intensity_column), which(is_intensity_column))
  x[, col_idx_new]
}


#' Update attribute 'intensity_columns' according to 'sam'
#'
#' @param x MS-DIAL alignment results table
#' @param sam sample table. If NULL, uses attr(x, "msdial_sam")
#'
#' @returns an object of the same class as 'x'
#' @noRd
updateIntensityColumnIndex <- function(x, sam = NULL) {
  attr_msdial_sam <- "msdial_sam"
  attr_intensity_columns <- "intensity_columns"
  sam_sample_name_column <- "name"
  if(is.null(sam))
    sam <- attr(x, attr_msdial_sam)
  col_idx <- match(sam[[ sam_sample_name_column ]], colnames(x))
  if(any(!is.finite(col_idx)))
    warning("'sam' not synchronized with column names")
  attr(x, attr_intensity_columns) <- col_idx
  x
}
