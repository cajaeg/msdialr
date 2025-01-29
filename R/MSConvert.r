#' Convert MS files to mzML, mzXML etc.
#'
#' Wrapper for ProteoWizard MSConvert command-line tool ("msconvert.exe").
#'
#' Supported input files are the same as in MSConvert GUI:
#' \itemize{
#' \item AB/Sciex T2D
#' \item Agilent MassHunter
#' \item Bruker BAF
#' \item mzML
#' \item MZ5
#' \item mzML
#' \item mzXML
#' \item Sciex WIFF(2)
#' \item Shimadzu LCD
#' \item Thermo RAW
#' \item Waters RAW
#' }
#'
#' Add \code{options(msdialr_msconvert_path = "/path/to/msconvert.exe")}
#' (adjusted accordingly to your system) to your .Rprofile file. Use
#' \code{Sys.getenv("R_USER")} to check where this file is located on your
#' system.
#'
#' Run \code{msconvert()} without arguments to see all conversion filter options
#' provided by MSConvert.
#'
#' @param infile file(s) to convert - see Details for supported file formats. If
#'   missing, show MSConvert help page.
#' @param outdir output directory or NULL for same as input file (default)
#' @param export_mode one of "mzML" (default), "mzXML" , "mzData", or "JCAMP"
#' @param binary_encoding numeric: 32 (default) or 64
#' @param filter character vector of filters to apply, "" for no filter
#'   (default)
#' @param n_cores number of CPU cores to use. If >1, the function will run in
#'   parallel mode using \code{parallel::makeCluster()}.
#'
#' @return character vector listing names of converted files
#' @export
#' @references Chambers, Matthew C, Brendan Maclean, Robert Burke, Dario Amodei,
#'   Daniel L Ruderman, Steffen Neumann, Laurent Gatto, et al. ‘A Cross-Platform
#'   Toolkit for Mass Spectrometry and Proteomics’. Nature Biotechnology 30, no.
#'   10 (October 2012): 918–20. https://doi.org/10.1038/nbt.2377.
#'
#' @examples
#'
#' # example for converting file
#' \dontrun{
#' mzml_file <- list.files(
#'   system.file("extdata", package = "msdialr"),
#'   "GCQTOF_PAH.*\\.mzML",
#'   full.names = TRUE
#' )
#' workdir <- tempfile() # you would usually convert in place, here we need to copy to temporary folder
#' dir.create(workdir)
#' file.copy(mzml_file, workdir)
#' new_path <- file.path(workdir, basename(mzml_file))
#' converted <- msconvert(new_path,
#'                        export_mode = "mzXML",
#'                        filter = "threshold absolute 10000 most-intense")
#' file.info(converted)
#'
#' unlink(workdir, recursive = TRUE, force = TRUE) # remove temporary folder
#' }
#'
msconvert <- function(infile,
                      outdir = NULL,
                      export_mode = c("mzML", "mzXML", "mz5", "mgf")[1],
                      binary_encoding = c(32, 64)[1],
                      filter = c("", "peakPicking true 1-", "threshold absolute 100 most-intense")[1],
                      n_cores = 1) {
  format_output <- function(x) {
    sub("writing output file: ", "", x[grep("^writing output file", x)])
  }
  if (missing(infile)) {
    do_msconvert() # display help
    rv <- 0
  } else {
    stopifnot(all(file.exists(infile)))
    n_files <- length(infile)
    rv <-
      if (n_files > 1) {
        if (n_cores == 1) {
          unlist(
            lapply(
              infile,
              msconvert,
              outdir = outdir,
              export_mode = export_mode,
              binary_encoding = binary_encoding,
              filter = filter,
              n_cores = 1
            )
          )
        } else {
          cl <- parallel::makeCluster(min(n_cores, n_files))
          on.exit(parallel::stopCluster(cl))
          unlist(
            parallel::parLapply(
              cl,
              infile,
              msdialr::msconvert,
              outdir = outdir,
              export_mode = export_mode,
              binary_encoding = binary_encoding,
              filter = filter,
              n_cores = 1
            )
          )
        }
      } else {
        infile_arg <- shQuote(normalizePath(infile))
        if (is.null(outdir)) {
          outdir <- dirname(infile)
        }
        outdir_arg <- c("--outdir", shQuote(normalizePath(outdir)))
        outfile_arg <- c("--outfile", shQuote(paste0(
          sub("\\.[^\\.]+$", "", basename(infile)), ".", export_mode
        )))
        export_mode_arg <- paste0("--", switch(
          export_mode,
          mzXML = "mzXML",
          mzML = "mzML",
          mz5 = "mz5",
          mgf = "mgf",
          ... = "mzML"
        ))
        binary_encoding_arg <- paste0("--", binary_encoding)
        filter_arg <- if (filter[1] != "") {
          unlist(lapply(filter, function(x)
            c("--filter", shQuote(x))))
        } else
          ""
        final_args <- c(
          infile_arg,
          outdir_arg,
          outfile_arg,
          export_mode_arg,
          binary_encoding_arg,
          filter_arg
        )
        format_output(do_msconvert(args = final_args))
      }
  }
  return(rv)
}


do_msconvert <- function(args = "--help") {
  path <- normalizePath(getOption("msdialr_msconvert_path"))
  if (!file.exists(path)) {
    stop("MSConvert not found, see ?msconvert for how to set path")
  }
  system2(path,
          args = args,
          stdout = TRUE,
          stderr = FALSE)
}


#' List raw and converted MS files
#'
#' @param path file path
#' @param suffix file suffix of original MS files
#' @param target_format file format of converted files
#'
#' @return data.frame
#' @export
checkMSfiles <- function(path,
                         suffix = c(".D", ".wiff"),
                         target_format = "mzML") {
  stopifnot(dir.exists(path))
  target_format <- match.arg(target_format, c("mzML"))
  ms_files <- lapply(
    sprintf("\\%s$", suffix),
    list.files,
    path = path,
    full.names = TRUE,
    ignore.case = TRUE
  )
  suffix_sel <- which.max(lengths(ms_files))
  ms_files <- ms_files[[suffix_sel]]
  target_files <- file.path(path, sprintf("%s.%s", sub(
    sprintf("\\%s$", suffix[suffix_sel]), "", basename(ms_files)
  ), target_format))
  data.frame(ms_file = ms_files,
             target_file = target_files,
             is_converted = file.exists(target_files))
}