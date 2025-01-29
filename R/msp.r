#' Read MSP library files
#'
#' @param mspfile file name
#' @param GMDsynonfix Fix for GMD MSP file (remove "Synon: " prefix)
#' @param .progress Show a progress bar, see \code{\link[plyr]{ldply}()}.
#'
#' @return data.frame
#' @export
#'
#' @examples
#' msp_file <- system.file("extdata/MSP_demo.msp", package = "msdialr")
#' msp_lib <- readMSP(msp_file)
#' head(msp_lib)
readMSP <-
  function(mspfile,
           GMDsynonfix = FALSE,
           .progress = c("text", "none")[1])
  {
    if (!file.exists(mspfile))
      stop("Input file not found")
    txt <- readLines(mspfile)
    sec_start <- grep("^NAME: ", txt, ignore.case = T)
    sec_split <- rep(seq_along(sec_start), diff(c(sec_start, length(txt) +
                                                    1)))
    tmp <- split(txt, sec_split)
    .progress <- if (length(tmp) <= 50)
      "none"
    else
      .progress # disable PB for short MSP's
    plyr::ldply(tmp, function(x) {
      ## parse sections
      flt <- nchar(x) > 0 # remove blank lines
      x <- x[flt]
      is_metadata <- 1:length(x) <= grep("^num peaks:", x, ignore.case = T) # it's safe to assume that "Num Peaks:" comes last
      if (GMDsynonfix) {
        x[is_metadata] <- sub("^Synon: ", "", x[is_metadata])
      }
      metadata <- strsplit(x[is_metadata], ": ")
      df <- data.frame(
        matrix(
          sapply(metadata, "[", 2),
          nrow = 1,
          ncol = sum(is_metadata),
          dimnames = list(NULL, c(sapply(metadata, "[", 1)))
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      convertsToNumeric <- which(!apply(df, 2, function(x)
        methods::is(tryCatch(
          as.numeric(x),
          warning = function(w)
            w
        ), "warning")))
      for (i in convertsToNumeric)
        df[, i] <- as.numeric(df[, i])
      pks <- x[!is_metadata]
      pks_str <- paste(gsub("\\t", ":", x[!is_metadata]), collapse = " ")
      df$PEAKS <- pks_str
      df
    }, .id = NULL, .progress = .progress)
  }


#' Write MSP library files
#'
#' @param x data.frame to be exported as MSP
#' @param file MSP file to be written
#' @param append append to existing file?
#'
#' @return None (invisible NULL)
#' @export
#'
#' @examples
#' msp <- data.frame(
#'   name = "Test_compound",
#'   num_peaks = 3,
#'   peaks = "100:999 101:122 102:33"
#' )
#' msp
#' writeMSP(msp, file = "")
#' \dontrun{
#' writeMSP(msp, file = "test.msp")
#' }
writeMSP <- function(x, file = NULL, append = FALSE) {
  stopifnot(inherits(x, "data.frame"))
  peak_columns <- lapply(c("^num[ _]peaks$", "^peaks$"), grep, tolower(colnames(x)))
  if (length(peak_columns) != 2 ||
      any(lengths(peak_columns)) != 1)
    stop("peak columns should be named 'Num Peaks' and 'PEAKS'")
  num_peak_column <- peak_columns[[1]]
  peak_columns <- peak_columns[[2]]
  meta_columns <- (1:ncol(x))[-peak_columns]
  nr <- nrow(x)
  cn.out <- toupper(colnames(x))
  cn.out[num_peak_column] <- gsub("_", " ", cn.out[num_peak_column])
  cidx1 <- meta_columns
  v1 <- paste(cn.out[cidx1], unlist(t(x[, cidx1])), sep = ": ")
  v1 <- gsub(": NA$", ":", v1)
  v1 <- split(v1, rep(1:nr, each = length(cidx1)))
  cidx2 <- peak_columns
  pks <- trimws(
    gsub(
      "[^0-9\\.]*([0-9\\.]+)[^0-9\\.]+([0-9\\.]+)[^0-9\\.]*",
      "\\1:\\2 ",
      trimws(x[[cidx2]])
    )
  )
  pks <- strsplit(x[[cidx2]], " ")
  v2 <- lapply(pks, sub, pattern = ":", replacement = " ")
  stopifnot(length(v1) == length(v2))
  v3 <- rep(list(""), length(v1))
  v <- mapply(c, v1, v2, v3, SIMPLIFY = FALSE)
  cat(
    unlist(v),
    file = if (!is.null(file))
      file
    else
      "",
    sep = "\n",
    append = append
  )
}
