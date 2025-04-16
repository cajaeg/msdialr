#' Annotate mass spectra using NIST MSPepSearch
#'
#' Use MSPepSearch command line tool to search spectral libraries. Download and
#' unpack the MSPepSearch zip file from NIST site, copy all files to "MSSEARCH"
#' directory and add the following line to
#' your .Rprofile: \code{options(msdialr_mspepsearch_path =
#' "/path/to/MSPepSearch64.exe")}
#'
#' Tested with MSPepSearch64 v2019_02_22.
#' @param spec mass spectrum: matrix or dataframe with m/z and intensity in the
#'   first two columns, or a list of matrices
#' @param searchtype one of the \emph{low-res} search types:
#' \describe{
#' \item{I}{identity, default}
#' \item{Q}{quick identity}
#' \item{S}{simple similarity}
#' \item{H}{hybrid similarity}
#' \item{L}{neutral loss similarity}
#' \item{M}{ms/ms in EI library}
#' }
#' ... or a \emph{high-res} search types:
#' \describe{
#' \item{P}{peptide ms/ms search in a ms/ms library (use peak annotations and weighting}
#' \item{G}{generic ms/ms search in a ms/ms library (no peak annotations or weighting used}
#' \item{D}{score is pure dot product}
#' }
#' @param presearchtype one of the following:
#' \describe{
#' \item{d}{standard (use pre-indexed peaks)}
#' \item{s}{sequential (compare to all library spectra)}
#' \item{f}{fast presearch (use pre-indexed peaks)}
#' \item{k\[n\]}{Compare spectra with the same n first segments of InChIKey (n=1,2,3); k is the same as k1}
#' \item{m}{precursor ion m/z within precursor ion m/z uncertainty - only for HiRes search \{PGD\}, incompatible with u, y}
#' }
#' @param searchoption one or more of the following as one character string, e.g. "rpx":
#' \describe{
#' \item{r}{reverse search}
#' \item{p}{penalize rare compounds (I in Mainlib & Replib only)}
#' \item{x}{do not output hit probabilities}
#' }
#' for high-res search (PGD) only:
#' \describe{
#' \item{j}{hit reJection; Pj: reject non-peptide hits, Gj: reject peptide hits}
#' \item{a}{alternative peak matching in computing match factors (recommended)}
#' \item{v}{output Rev-Dot (reverse match factor)}
#' \item{i}{ignore peaks within zi or max(z,m) m/z units around precursor m/z; see /ZI, /Z, /M option}
#' \item{leh\[n\]}{ms/ms search threshold: l Low, e mEdium, h High(default), hn (n=1..9) Higher}
#' }
#' for peptide search (P) only:
#' \describe{
#' \item{oqn}{o Use OMSSA score; q Q-tof; n use number of replicates}
#' }
#' @param libs libraries to search (default: c("mainlib", "replib"))
#' @param hits how many hits to return, default 10
#' @param precursor_mz_dev /Z\[PPM\] parameter = precursor ion m/z uncertainty,
#'   0.00006-500 m/z units (default=1.6) or 0.015-1e5 ppm, only for {PGD}
#' @param ms_mz_dev /M\[PPM\] parameter = ms peak m/z uncertainty, 0.00006-500 m/z units (default=0.6) or 0.015-1e5 ppm, only for {PGD}
#' @param mzlimits set m/z limits for spectral comparison
#' @param ri_obs measured retention index (RI)
#' @param max_ri_dev passed to \code{scoreRI()}
#' @param ri_minscore minimum RI score; passed to \code{scoreRI()}

#' @param localRdata attach meta data from local .rda (path to .rda or object
#'   already loaded in workspace)
#' @param .debug if TRUE keep temporary files created by MSPepSearch and return raw MSPepSearch output
#' @return data.frame
#' @export
#'
#' @examples
#' # standard nominal mass EI search
#' \dontrun{
#' spec <- str2spec("39:74 41:150 65:49 77:84 91:72 95:142 107:275 135:999 136:103 150:262")
#' mspepsearch(spec)
#' mspepsearch(spec, ri_obs = 1300)
#' }
#' # high-res search without and with precursor m/z
#' \dontrun{
#' spec <- str2spec(
#'   "57.0696:128 69.0332:101 86.0597:146 97.0282:641 114.0548:307 125.0595:86 129.0657:197 157.0605:999 185.0918:168 213.1231:76"
#' )
#' mspepsearch(spec,
#'             searchtype = "G",
#'             searchoption = "u",
#'             libs = "nist_msms")
#' mspepsearch(
#'   spec,
#'   searchtype = "G",
#'   searchoption = "m",
#'   precursor_mz = 213.1234,
#'   libs = "nist_msms"
#' )
#' }
#' # peptide search (= weighted high-res search)
#' \dontrun{
#' spec <- str2spec(
#'   "226.0814:179 243.108:999 342.1763:539 411.1969:194 519.2925:326 524.2802:225 620.3403:345 820.4613:225 1161.6263:473 1162.6259:263"
#' )
#' mspepsearch(
#'   spec,
#'   searchtype = "P",
#'   searchoption = "m",
#'   precursor_mz = 1161.63,
#'   libs = "nistmab_v20190711"
#' )
#' }
mspepsearch <- function(spec,
                        searchtype = "I",
                        presearchtype = "d",
                        searchoption = "v",
                        libs = c("mainlib", "replib"),
                        hits = 10,
                        precursor_mz = NULL,
                        precursor_mz_dev = NULL,
                        ms_mz_dev = NULL,
                        mzlimits = c(-1, -1),
                        ri_obs = NULL,
                        max_ri_dev = 100,
                        ri_minscore = 0.9,
                        localRdata = NULL,
                        .debug = FALSE) {
  ri_dev = 5 # argument removed from argument list as it doesn't affect MSPepSearch scores; use max_ri_dev instead
  # previous RD text: allowed RI tolerance (1-32767)
  ri_penalty = 5 # argument removed from argument list as it doesn't have affect MSPepSearch scores; use ri_minscore instead
  # previous RD text: penalty rate (1-65535); NIST presets: Very weak(10), WK=Weak(20), AV=Average(50), ST=Strong(100), VS=Very strong(200), IN=Infinite(65535)
  mspepsearch <- getOption("msdialr_mspepsearch_path")
  if (is.null(mspepsearch) || !file.exists(mspepsearch)) {
    stop("Option \"msdialr_mspepsearch_path\" not set or tool not found. See help file.")
  }
  nistmssearchdir <- dirname(mspepsearch)
  format_cas <- function(CAS) {
    CAS <- as.character(CAS)
    if (CAS == "0")
      return(NA_character_)
    sprintf(
      "%s-%s-%s",
      substr(CAS, 1, nchar(CAS) - 3),
      substr(CAS, nchar(CAS) - 2, nchar(CAS) - 1),
      substr(CAS, nchar(CAS), nchar(CAS))
    )
  }
  ## --- check arguments ---
  if (missing(spec)) {
    message("No spectrum provided, displaying MSPepSearch help page.")
    owd <- setwd(dirname(mspepsearch))
    on.exit(setwd(owd), add = TRUE)
    system2(command = basename(mspepsearch))
    return(invisible(NULL))
  }
  if (is.data.frame(spec) || is.matrix(spec)) {
    spec <- list(spec)
  }
  if (!is.null(ri_obs) && length(ri_obs) != length(spec))
    ri_obs <- rep(ri_obs[1], length(spec))
  stopifnot(all(libs %in% tolower(basename(
    list.dirs(nistmssearchdir, recursive = FALSE)
  ))))
  stopifnot(
    length(searchtype) == 1 &&
      nchar(searchtype) == 1 &&
      searchtype %in% c("I", "Q", "S", "H", "L", "M", "P", "G", "D")
  )
  stopifnot(
    length(presearchtype) == 1 &&
      substr(presearchtype, 1, 1) %in% c("d", "f", "k", "m", "s")
  )
  stopifnot(is.character(searchoption))
  if(presearchtype == "m" && is.null(precursor_mz)) {
    stop("precursor_mz required for presearchtype == \"m\"")
  }
  ## --- set tempdir and txt files ---
  tmpdir <- tempdir()
  tmpdir_work <- file.path(tmpdir, "searchNIST")
  if (!.debug) {
    on.exit(unlink(tmpdir_work, recursive = TRUE, force = TRUE))
  }
  dir.create(tmpdir_work, recursive = TRUE)
  mspfile <- file.path(tmpdir_work, "tmp.msp")
  outfile <- file.path(tmpdir_work, "tmp.txt")
  ##
  ## --- write msp file ---
  ##
  if (file.exists(mspfile))
    unlink(mspfile, recursive = F, force = T)
  n_spec <- length(spec)
  mspNames <- sprintf("RQuery_%06d", 1:n_spec)
  showProgress <- n_spec > 100
  if (showProgress) {
    message("Exporting spectra to msp file ...")
    pb <- utils::txtProgressBar(min = 0,
                                max = length(mspNames),
                                style = 3)
  }
  zz_len <- 1e6 # size of pre-allocated character vector
  zz <- vector("character", zz_len)
  j = 1
  for (i in 1:length(spec)) {
    zz[j] <- sprintf("Name: %s", mspNames[i])
    j = j + 1
    if (!is.null(ri_obs)) {
      zz[j] <- sprintf("Retention_index: SemiStdNP=%d/1/4 StdNP=%d",
                       round(ri_obs[i]),
                       round(ri_obs[i]))
      j = j + 1
    }
    if (!is.null(precursor_mz)) {
      zz[j] <- sprintf("PrecursorMZ: %.4f", precursor_mz[i])
      j = j + 1
    }
    m <- as.matrix(spec[[i]][, 1:2])
    zz[j] <- sprintf("Num Peaks: %d", nrow(m))
    j = j + 1
    tmp <- trimws(apply(m, 1, paste, collapse = " "))
    for (k in 1:length(tmp)) {
      zz[j] <- tmp[k]
      j = j + 1
    }
    zz[j] <- ""
    j = j + 1
    if (j > zz_len * 0.9) {
      cat(zz[1:j],
          file = mspfile,
          sep = "\n",
          append = file.exists(mspfile)) # space gets tight, write batch to disk
      zz <- vector("character", zz_len)
      j = 1
    }
    if (showProgress)
      utils::setTxtProgressBar(pb, i)
  }
  if (showProgress)
    close(pb)
  cat(zz[1:j],
      file = mspfile,
      sep = "\n",
      append = file.exists(mspfile))
  ##
  ## --- create argument list ---
  ##
  args <- c(
    sprintf("%s%s%s", presearchtype, searchoption, searchtype),
    # "v" allows reverse MF output with /OUTTAB
    if ("mainlib" %in% libs)
      "/MAIN mainlib",
    if ("replib" %in% libs)
      "/REPL replib",
    if (!is.null(ri_obs))
      sprintf("/RI sut%dr%d", ri_dev, ri_penalty),
    if (any(!libs %in% c("mainlib", "replib")))
      paste(paste("/LIB", libs[!libs %in% c("mainlib", "replib")]), collapse = " "),
    sprintf("/INP %s", normalizePath(mspfile)),
    sprintf("/HITS %d", hits),
    sprintf("/MzLimits %d %d", mzlimits[1], mzlimits[2]),
    if (!is.null(precursor_mz_dev))
      sprintf("/Z %.4f", precursor_mz_dev),
    if (!is.null(ms_mz_dev))
      sprintf("/M %.4f", ms_mz_dev),
    sprintf("/COL mw,nn,nm,cf,cn,ik,pz,dz,tz,tq,ce"),
    sprintf("/OUTTAB %s", suppressWarnings(normalizePath(outfile)))
  )
  ##
  ## --- run command ---
  ##
  owd <- setwd(dirname(mspepsearch))
  on.exit(setwd(owd), add = TRUE)
  if (showProgress)
    message("Calling MSPepSearch ...")
  rv0 <- system2(
    command = basename(mspepsearch),
    args = args,
    stdout = TRUE,
    stderr = FALSE
  )
  rv <- iconv(rv0, "IBM860", "UTF-8") # fix character encoding
  if (showProgress)
    message("Reading MSPepSearch results ...")
  txt0 <- readLines(outfile)
  ## browser()
  ##
  ## --- post-process command output ---
  ##
  lidx <- grep("^>$", txt0) + c(2, -1)
  nr <- diff(lidx) + 1
  outcols <- c(
    "query",
    "rank",
    "name",
    "formula",
    "mf",
    "rmf",
    "prob",
    "score",
    "nummp",
    "charge",
    "precursor_mz",
    "precursor_type",
    "precursor_mz_dev",
    "ce",
    "mods",
    "pep",
    "protein",
    "tfqry",
    "ri_obs",
    "ri_lib",
    "ri_score",
    "total_score",
    "column_type",
    "match_no_ri",
    "cas",
    "mw",
    "exactmass",
    "lib",
    "lib_id",
    "nistrn",
    "inchikey",
    "peaks"
  )
  out <- data.frame(matrix(
    nrow = max(0, nr),
    ncol = length(outcols),
    dimnames = list(NULL, outcols)
  ),
  stringsAsFactors = FALSE)
  NISTcn <- strsplit(txt0[lidx[1] - 1], "\t")[[1]] # column names in NIST output
  hasRI <- length(grep("^RI$", NISTcn)) > 0
  if (nr > 0) {
    txt <- txt0[seq(lidx[1], lidx[2])]
    txt <- strsplit(txt, "\t")
    out$query <- as.numeric(sub("RQuery_([0-9]+).*", "\\1", sapply(txt, "[[", 1)))
    out$rank <- as.numeric(sapply(txt, "[[", grep("^Rank$", NISTcn)))
    out$formula <- trimws(sapply(txt, "[[", grep("^Formula$", NISTcn)))
    out$prob <- as.numeric(sapply(txt, "[[", grep("^Prob\\(", NISTcn)))
    out$cas <- sapply(txt, function(x) {
      format_cas(x[[grep("^CAS$", NISTcn)]])
    })
    out$mw <- as.numeric(sapply(txt, "[[", grep("^Lib MW$", NISTcn)))
    out$exactmass <- as.numeric(sapply(txt, "[[", grep("^Mass$", NISTcn)))
    out$lib <- trimws(sapply(txt, "[[", grep("^Library$", NISTcn)))
    out$lib_id <- as.numeric(sapply(txt, "[[", grep("^Id$", NISTcn)))
    out$nistrn <- as.numeric(sapply(txt, "[[", grep("^NIST r.n.$", NISTcn)))
    out$inchikey <- trimws(sapply(txt, "[[", grep("^InChIKey$", NISTcn)))
    if (searchtype %in% c("I", "Q", "S", "H", "L", "M")) { # --- nominal mass search types
      out$name <- trimws(sapply(txt, "[[", grep("^Name$", NISTcn)))
      out$mf <- as.numeric(sapply(txt, "[[", grep("^MF$", NISTcn)))
      out$rmf <- as.numeric(sapply(txt, "[[", grep("^R.Match$", NISTcn)))
    } else if (searchtype %in% c("P", "G", "D")) { # --- high-res search types
      out$name <- trimws(sapply(txt, "[[", grep("^Peptide$", NISTcn)))
      out$mf <- as.numeric(sapply(txt, "[[", grep("^Dot Product$", NISTcn)))
      out$score <- as.numeric(sapply(txt, "[[", grep("^Score$", NISTcn)))
      out$nummp <- as.numeric(sapply(txt, "[[", grep("^NumMP$", NISTcn)))
      out$charge <- as.numeric(sapply(txt, "[[", grep("^Charge$", NISTcn)))
      try({
        out$precursor_mz <- as.numeric(sapply(txt, "[[", grep("^Lib Precursor m/z$", NISTcn)))
      }, silent = TRUE)
      try({
        out$precursor_type <- trimws(sapply(txt, "[[", grep("^Prec\\.Type$", NISTcn)))
      }, silent = TRUE)
      if ("m" %in% searchoption) {
        try({
          out$precursor_mz_dev <- as.numeric(sapply(txt, "[[", grep("^Delta\\(m/z\\)$", NISTcn)))
        }, silent = TRUE)
      }
      try({
        out$mods <- trimws(sapply(txt, "[[", grep("^Mods$", NISTcn)))
      }, silent = TRUE)
      try({
        out$pep <- trimws(sapply(txt, "[[", grep("^Pep$", NISTcn)))
      }, silent = TRUE)
      try({
        out$protein <- trimws(sapply(txt, "[[", grep("^Protein$", NISTcn)))
      }, silent = TRUE)
      try({
        out$tfqry <- trimws(sapply(txt, "[[", grep("^T/F-qry$", NISTcn)))
      }, silent = TRUE)
      try({
        out$ce <- round(as.numeric(sub("eV", "", trimws(
          sapply(txt, "[[", grep("^CE$", NISTcn))
        ))))
      }, silent = TRUE
      )
    }
    if (!is.null(ri_obs))
      out$ri_obs <- round(rep(ri_obs, tapply(out$query, out$query, length)))
    if (hasRI) {
      try({
        out$ri_lib <- as.numeric(gsub("[^0-9]", "", sapply(txt, "[[", grep(
          "^RI$", NISTcn
        ))))
        out$ri_score <- round(
          scoreRI(
            out$ri_obs - out$ri_lib,
            max_ri_dev = max_ri_dev,
            minscore = ri_minscore
          ),
          2
        )
        out$total_score <- round(out$mf * out$ri_score)
        out$column_type <- sub(".*([A-Z])", "\\1", sapply(txt, "[[", grep("^RI$", NISTcn)))
        out$match_no_ri <- as.numeric(gsub("[^0-9]", "", sapply(
          txt, "[[", grep("^Match no RI$", NISTcn)
        )))
      }, silent = TRUE)
    }
  }
  out <- out[order(out$query, if (hasRI)
    - out$total_score
    else-out$mf), ]
  out$rank <- as.numeric(unlist(tapply(out$query, out$query, seq_along)))
  if (is.null(localRdata) && "mainlib" %in% libs) {
    localRdata <- getOption("msdialr_localNISTrda")
  }
  if (!is.null(localRdata)) {
    if (inherits(localRdata, "data.frame")) {
      nistdb <- localRdata
    } else if (file.exists(localRdata)) {
      load(localRdata) # "nistdb"
    }
    if (exists("nistdb")) {
      idx <- match(out$nistrn, nistdb$nistrn)
      if (sum(is.finite(idx)) > 0) {
        out$peaks[which(is.finite(idx))] <- nistdb$peaks[idx[is.finite(idx)]]
      }
    }
  }
  attr(out, "timeinfo") <- txt0[length(txt0)]
  if (.debug) {
    print(tmpdir_work)
    return(txt0)
  } else {
    return(out)
  }
}


#' Annotate mass spectra using NIST MSPepSearch
#'
#' @param x MS-DIAL data.frame as returned by \code{\link{loadConsoleResults}()} or \code{\link{loadAlignmentResults}()}
#' @param in_column name of column containing prepared spectra (see \code{\link{prepareSpectra}()})
#' @param ri_column name of column containing RI values or 'NULL' to exclude RI
#' @param out_column name of column added to data.frame
#' @param ... further arguments passed to \code{\link{mspepsearch}()}
#' @return 'x' with new column 'out_column'
#' @export
#' @examples
#' # see ?loadAlignmentResults()
searchNIST <- function(x,
                       in_column = "s",
                       ri_column = c("average_ri", "ri"),
                       out_column = "NISTres",
                       ...) {
  stopifnot(inherits(x, c("matrix", "data.frame", "tbl", "tbl_df")))
  stopifnot(in_column %in% colnames(x))
  spectra <- x[[ in_column ]]
  if(!is.null(ri_column) && any(ri_column %in% colnames(x))) {
    ri_obs <- x[[ ri_column[which(ri_column %in% colnames(x))[1]] ]]
  } else {
    ri_obs <- NULL
  }
  tmp <- mspepsearch(spec = spectra, ri_obs = ri_obs, ...)
  tmp$query <- factor(tmp$query, levels = 1:length(spectra))
  x[[ out_column ]] <- split(tmp, tmp$query)
  x |> 
    relocateIntensityColumns()
}


#' Copy NIST matches to main data table 
#'
#' @param x MS-DIAL data.frame as returned by \code{\link{loadConsoleResults}()}
#' @param in_column name of column containing NIST results (default: "NISTres")
#' @param minscore copy NIST match only when score is >= minscore
#' @param eval which score to evaluate, typically "total_score" when RI is used and "mf" otherwise
#'
#' @return same as 'x'
#' @export
#'
#' @examples
#' # see ?loadAlignmentResults()
acceptNISTres <- function(x, 
                          in_column = "NISTres",
                          minscore = 700, 
                          eval = c("total_score", "mf", "rmf")) {
  out <- 
    if(in_column %in% colnames(x)) {
      for(i in 1:nrow(x)) {
        nr <- x[[ in_column ]][[ i ]]
        eval_col <- eval[sapply(eval, function(x) any(is.finite(nr[[x]])))][1] # use first valid option in 'eval'
        nr_flt <- nr[nr[[ eval_col ]] >= minscore, ]
        if(nrow(nr_flt) > 0) {
          nr_flt <- nr_flt[ order(nr_flt[[ eval_col ]], decreasing = TRUE), ]
          x$metabolite_name[i] <- nr_flt$name[1]
          x$formula[i] <- nr_flt$formula[1]
          x$reference_ri[i] <- nr_flt$ri_lib[1]
          x$ri_similarity[i] <- nr_flt$ri_score[1]
          x$inchikey[i] <- nr_flt$inchikey[1]
          x$dot_product[i] <- nr_flt$mf[1]
          x$reverse_dot_product[i] <- nr_flt$rmf[1]
          x$total_score[i] <- nr_flt$total_score[1]
        }
      }
      x
    } else {
      message("no NIST results found - returning input unchanged")
      x
    }
  out |>
    relocateIntensityColumns() 
}


#' Calculate a simple score based on RI deviation
#'
#' @param ri_dev deviation of observed RI from library RI
#' @param max_ri_dev deviation greater than this receive the min score; NA values also receive the min score
#' @param minscore minimum score, default 0.75
#' @param maxscore maximum score, default 1.0
#' @return numeric vector
#' @export
#'
#' @examples
#' ri_obs <- 1225
#' ri_lib <- 1200
#' scoreRI(ri_obs - ri_lib)
scoreRI <- function(ri_dev,
                    max_ri_dev = 100,
                    minscore = 0.75,
                    maxscore = 1) {
  ri_dev[!is.finite(ri_dev)] <- max_ri_dev
  ri_dev <- abs(ri_dev)
  ((pmin(ri_dev, max_ri_dev) - 0) / (max_ri_dev - 0)) * (minscore - maxscore) + maxscore
}

