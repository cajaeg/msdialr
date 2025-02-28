#' check if arg is present in arg list, if not, set it to 'val' 
#'
#' (internal function)
#' @param args arg list
#' @param name arg
#' @param val default value
#'
#' @return list
#' @noRd
#'
checkArgs <- function(args, name, val) {
  if(!name %in% names(args))
    args[[name]] <- val
  return(args)
}


#' Check if a package is installed
#'
#' (internal function)
#' @param pkg package to check (character(1))
#' @param libPaths library paths to check, defaults to \code{.libPaths()}
#'
#' @return logical
#' @noRd
#' 
isInstalled <- function(pkg, libPaths = .libPaths()) {
  libPaths <- match.arg(libPaths)
  any(grepl(pkg, basename(
    list.dirs(libPaths, recursive = FALSE)
  )))
}

#' Numeric equivalent to 'match'
#' 
#' @param x numeric vector
#' @param table numeric vector against which to match x
#' @param delta_x allowed absolute difference
#' @return array index
#' @export
#'
#' @examples
#' x <- c(3:5)+0.01
#' y <- 1:10
#' nummatch(x, y, delta_x = 0.1)
nummatch <- function(x, table, delta_x = 0.005) {
  do_match <- function(x, table, delta_x) {
    which(abs(table - x) < delta_x)
  }
  pidx <- mapply(do_match, x, MoreArgs = list(table=table, delta_x=delta_x))
  return(cbind(
    x = rep(seq_along(pidx), sapply(pidx, length)),
    table = unlist(pidx)
  ))
}

#' m/z or RT clustering based on 'hclust'
#'
#' @param x numeric vector
#' @param h tree height at which to cut, see \code{\link[stats]{cutree}()}
#' @param k desired number of groups, see \code{\link[stats]{cutree}()}
#' @param method cluster method, see \code{\link[stats]{hclust}()}
#' @param useFastCluster use \code{hclust.vector()} from \code{fastCluster}
#'   package instead of the standard \code{stats::hclust()}. Default is to use
#'   \code{fastCluster} if it is installed.
#'
#' @return numeric vector representing group structure
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- abs(rnorm(100, 30, 15))
#' g <- hcgroup(x, h = 5, useFastCluster = FALSE)
#' g1 <- hcgroup(x, h = 5, useFastCluster = TRUE)
#' identical(g, g1)
#' plot(data.frame(x = x, y = 1), type = "h", col = g, ylim = c(0, 1))
hcgroup <- function(x,
                    h = NULL,
                    k = NULL,
                    method = "centroid",
                    useFastCluster = NULL) {
  smallNumber <- .Machine$double.xmin
  if (is.null(useFastCluster) ||
      (is.logical(useFastCluster) && useFastCluster))
    useFastCluster <- isInstalled("fastcluster")
  if (length(x) >= 2) {
    cl <- if (useFastCluster) {
      fastcluster::hclust.vector(c(smallNumber, x), method = method)
    } else {
      hc <- stats::hclust(stats::dist(c(smallNumber, x)) ^ 2, method = method)
      hc$height <- sqrt(hc$height)
      hc
    }
    tryCatch(
      stats::cutree(cl, h = h, k = k)[-1] - 1, 
      error = function(e) {
        cl$height <- round(cl$height, 6)
        stats::cutree(cl, h = h, k = k)[-1] - 1
      })
  }
  else
    1
}


#' replace zero values with NA
#'
#' (internal function)
#' @param x numeric vector
#'
#' @return numeric vector
#' @noRd
#'
zero2NA <- function(x) {
  x[x==0] <- NA
  return(x)
}


#' replace NA with zero
#'
#' (internal function)
#' @param x numeric vector
#'
#' @return numeric vector
#' @noRd
#'
NA2zero <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}


#' replace Inf with zero
#'
#' (internal function)
#' @param x numeric vector
#'
#' @return numeric vector
#' @noRd
#'
Inf2zero <- function(x) {
  x[!is.finite(x)] <- 0
  return(x)
}


#' simple replacement for 'plyr::ldply()'
#'
#' (internal function)
#' @param l list
#'
#' @return data.frame
#' @noRd
#'
do.rbind <- function(l) {
  data.frame(do.call(rbind, l), id = rep(1:length(l), sapply(l, nrow)))
}


#' Skip plot in multi-figure layout by plotting nothing
#'
#' @param message optional text message to show, e.g. "no data"
#'
#' @return NULL
#' @export
#'
emptyplot <- function(message = NULL) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 2), ylim = c(0, 2))
  if(!is.null(message))
    graphics::text(1, 1, labels = message)
  invisible(NULL)
}


#' Calculate m/z from chemical formula and adduct
#'
#' @param x chemical formula
#' @param adduct adduct ("\[M+H\]+", "\[M\]+" etc.) or NULL for neutral mass
#' @param round digits to round m/z to
#'
#' @return data.frame
#' @export
#'
#' @examples
#' fml2mz("C6H6O", "[M]+") # (radical) cation
#' fml2mz("C6H6O") # neutral mass
#' fml2mz("C3H7NO2", c("[M+H]+", "[M-H]-")) # (de)protonated ion
fml2mz <- function(x, adduct = NULL, round = 4) {
  isotopes <- NULL
  utils::data("isotopes", package = "enviPat", envir = environment()) # isotopes
  if (is.character(x)) {
    fml <- x
    m <- base::round(enviPat::check_chemform(isotopes, fml)[, "monoisotopic_mass"], round)
  } else if (is.numeric(x)) {
    fml <- NA_character_
    m <- x
  } else {
    return(NA)
  }
  if (is.null(adduct)) {
    r <- data.frame(
      name = "M",
      nmol = 1,
      charge = 0,
      massdiff = 0
    )
  } else {
    r <- ion2rule(adduct)
  }
  adductmz <- (m * r$nmol + r$massdiff) / pmax(1, abs(r$charge))
  n_adducts <- nrow(r)
  data.frame(
    formula = rep(fml, n_adducts),
    neutral_mass = rep(m, n_adducts),
    adduct = r$name,
    nmol = r$nmol,
    charge = r$charge,
    mz = base::round(adductmz, round)
  )
}


#' Calculate isotope pattern using 'enviPat'
#' 
#' (internal function)
#' @param fml chemical formula
#' @param resolution mass spectral resolution (m/delta_m)
#' @param ... passed to 'enviPat::isowrap()'
#'
#' @return list
#' @noRd
#'
#' @examples
#' \dontrun{
#' fml2iso("C6H6O")
#' }
fml2iso <- function(fml, resolution = 20000, ...) {
  isotopes <- NULL
  utils::data("isotopes", package = "enviPat", envir = environment()) # "isotopes"
  fixMS <- function(ms) {
    colnames(ms) <- c("mz", "i")
    ms[, 1] <- round(ms[, 1], 5)
    ms[, 2] <- round(ms[, 2], 1)
    ms
  }
  out <- suppressMessages(
    enviPat::isowrap(
      isotopes,
      enviPat::check_chemform(isotopes, fml),
      resmass = FALSE,
      resolution = resolution,
      ...
    )
  )
  return(lapply(out, fixMS))
}


#' Derive mass difference from ion notation
#'
#' (internal function used by fml2mz())
#' @param ions character vector, e.g. "\[M+H\]+"
#'
#' @return data.frame
#' @noRd
#'
#' @examples
#' \dontrun{
#' ion2rule("[M+H]+")
#' }
ion2rule <- function(ions = "[M+H]+") {
  checkSymbol <- function(ion) {
    regexpr("\\[[0-9]{0,2}M.*\\][0-9]{0,2}[\\+\\-]{1,2}", ion) != -1
  }
  shortCuts <- cbind(
    c("M+H", "M+Na", "M+K", "M+NH4", "M+", "M", "M.", "M-H", "M+Cl-", "M-"),
    c(
      "[M+H]+",
      "[M+Na]+",
      "[M+K]+",
      "[M+NH4]+",
      "[M]+",
      "[M]+",
      "[M]+",
      "[M-H]-",
      "[M+Cl]-",
      "[M]-"
    )
  )
  emass <- NULL
  chemical_elements <- NULL
  utils::data(emass, envir = environment(), package = "msdialr")
  utils::data(chemical_elements, envir = environment(), package = "msdialr")
  out <- lapply(ions, function(ion) {
    if (ion %in% shortCuts[, 1])
      ion <- shortCuts[, 2][which(shortCuts[, 1] == ion)]
    if (!checkSymbol(ion))
      stop("invalid ion")
    nmol <- sub(".*[^0-9M]([0-9]?M).*", "\\1", ion)
    nmol <- sub("M", "", nmol)
    nmol <- as.numeric(ifelse(nmol == "", 1, nmol))
    ch <- sub(".*[^0-9]([0-9]{0,2}[\\+\\-])$", "\\1", ion)
    sgn <- sub("[^\\+\\-]", "", ch)
    sgn <- ifelse(sgn == "+", 1, -1)
    ch <- sub("[\\+\\-]", "", ch)
    ch[ch == ""] <- "1"
    ch <- as.numeric(ch)
    ch <- ch * sgn
    x <- ion
    x <- sub("^.*\\[", "", x)
    x <- sub("\\].*", "", x)
    x <- sub("[0-9]?M", "", x)
    starts <- gregexpr("[\\+\\-]", x)[[1]]
    ends <- c(starts[-1] - 1, nchar(x))
    n <- length(starts)
    spl <- lapply(1:n, function(i)
      substr(x, starts[i], ends[i]))
    massdiff <- lapply(spl, function(y) {
      sgn <- sub("^([\\+\\-]).*", "\\1", y)
      sgn <- ifelse(sgn == "+", 1, -1)
      el <- sub("^[\\+\\-]", "", y)
      if (regexpr("^[0-9]+[A-Za-z]+", el) != -1)
        el <- gsub("([0-9]+)([A-Za-z]+)", "\\2\\1", el)
      el <- fml2tbl(el)
      masses <- sapply(colnames(el), function(a) {
        chemical_elements[, 2][which(chemical_elements[, 1] == a)[1]]
      })
      return(sum(masses * el[1, ]) * sgn)
    })
    massdiff <- sum(unlist(massdiff), na.rm = TRUE) + ch *
      -emass
    return(
      data.frame(
        name = ion,
        nmol = nmol,
        charge = ch,
        massdiff = massdiff,
        stringsAsFactors = FALSE
      )
    )
  })
  return(do.call("rbind", out))
}


#' Tabulate elements in chemical formula
#'
#' (internal function used by fml2mz())
#' @param fml chemical formula
#' @param elements chemical elements to consider, or NULL for all elements
#'   occurring in 'fml'
#'
#' @return data.frame
#' @noRd
#' 
#' @examples
#' \dontrun{
#' fml2tbl("C6H6O")
#' fml2tbl("C6H6O", elements = c("C", "H", "O", "S"))
#' }
#' 
fml2tbl <- function (fml, elements = NULL) {
  getElements <- function(form) {
    starts <- gregexpr("[A-Z]", form)[[1]]
    stops <- c(starts[-1] - 1, nchar(form))
    elements <- sapply(1:length(starts), function(i)
      substr(form, starts[i], stops[i]))
    cnts <- sub("[A-Za-z]+", "", elements)
    cnts[cnts == ""] <- 1
    cnts <- as.numeric(cnts)
    elements <- sub("[0-9]+", "", elements)
    return(sort(rep(elements, cnts)))
  }
  checkFormula <- function(form) {
    el <- getElements(form)
    paste0(names(table(el)), table(el), collapse = "")
  }
  if (!is.character(fml))
    return(NULL)
  fml <- gsub("[^A-Za-z0-9]", "", fml)
  fml <- checkFormula(fml)
  if (is.null(elements)) {
    elements <- unique(unlist(lapply(fml, getElements)))
  }
  out <- matrix(ncol = length(elements), nrow = length(fml))
  for (i in seq(along = elements)) {
    el <- elements[i]
    rgx <- sprintf(".*%s([0-9]+).*", el)
    cnt <- sub(rgx, "\\1", fml)
    cnt[regexpr("[[:alpha:]]", cnt) == 1] <- "0"
    cnt[cnt == ""] <- "1"
    cnt <- as.numeric(cnt)
    out[, i] <- cnt
  }
  colnames(out) <- elements
  return(out[, order(colnames(out)), drop = FALSE])
}


#' Find last created file
#'
#' @param path file path
#' @param pattern file pattern
#' @param ... passed to \code{list.files()}
#'
#' @returns character()
#' @noRd
findLatestFile <- function(path = ".", pattern = NULL, ...) {
  fls <- list.files(path = path, pattern = pattern, full.names = TRUE, ...)
  if(length(fls) > 0) {
    mtime <- file.info(fls)$mtime
    sel <- order(mtime, decreasing = TRUE)[1]
    message(sprintf("Found %d file(s), selecting %s", length(fls), basename(fls[sel])))
    fls[sel]
  } else {
    fls
  }
}


#' Find last created alignment results file in a given folder
#'
#' @param path file path
#' @param pattern file pattern, defaults to "Height.*_\[0-9\]+\\.txt"
#' @param ... passed to \code{list.files()}
#'
#' @returns character()
#' @export
findLatestHeightFile <- function(path = ".", pattern = "Height.*_[0-9]+\\.txt", ...) {
  findLatestFile(path = path, pattern = pattern, ...)
}


#' Format numeric IDs
#'
#' @param x numeric
#' @param prefix prefix
#' @param n_digits n_digits
#'
#' @returns character()
#' @export
#' @examples
#' formatNumericIDs(1:10)
#' formatNumericIDs(1:1000)[1:10]
formatNumericIDs <- function(x,
                             prefix = "",
                             n_digits = NULL) {
  if(is.null(n_digits))
    n_digits <- ceiling(log10(max(x))) + 1
  sprintf(paste0(prefix, "%0", n_digits, "d"), x)
}


#' Convert m/z to neutral mass
#'
#' @param mz numeric()
#' @param adduct character()
#'
#' @returns numeric()
#' @noRd
mz2nm <- function(mz, adduct = "[M+H]+") {
  r <- ion2rule(adduct)
  nmol <- r[, "nmol"]
  z <- r[, "charge"]
  delta_mass <- r[, "massdiff"]
  (mz * abs(z) - delta_mass) / nmol
}


#' Convert mass spectrum from character to matrix representation
#'
#' @param x mass spectrum represented as character string, e.g. "100:999 101:121"
#' @param colnames column names in resulting data.frame (default \code{c("mz",
#'   "i")})
#'
#' @return data.frame
#' @export
#'
#' @examples
#' str2spec("100:999 101:121")
str2spec <- function (x, colnames = c("mz", "i"))
{
  tmp <- trimws(gsub("[^0-9\\.]*([0-9\\.]+)[^0-9\\.]+([0-9\\.]+)[^0-9\\.]*", "\\1 \\2 ", trimws(x)))
  tmp <- trimws(unlist(strsplit(tmp, " ")))
  tmp <- as.numeric(tmp)
  is_mz <- 1:length(tmp) %% 2 == 1
  out <- data.frame(mz = tmp[is_mz], i = tmp[!is_mz])
  base::colnames(out) <- colnames
  return(out)
}


#' Extract m/z portion of MS character string representation
#'
#' @param x character string
#' @param split character vector
#'
#' @return numeric vector of m/z values
#'
str2mz <- function(x, split = " ") {
  as.numeric(unlist(strsplit(x, split)[[1]]))
}


#' Convert mass spectrum from matrix to character representation
#'
#' @param x two-column matrix or data.frame
#' @param sep separator to be used between 'mz' and 'intensity' (default ':')
#' @param collapse separator to be used between (m/z, i) blocks (default ' ')
#' @param round Number of decimal places to be used for mz and intensity,
#'   respectively (default \code{c(4,0)})
#'
#' @return character string
#' @export
#'
#' @examples
#' s <- data.frame(mz = c(100, 101), i = c(999, 121))
#' spec2str(s)
spec2str <- function(x, sep = ":", collapse = " ", round = c(4,0)) {
  if(!any(inherits(x, c("matrix", "data.frame"), which = TRUE) == 1)) {
    stop("no matrix or data.frame provided")
  }
  if(!is.null(round)) {
    x[,1] <- base::round(x[,1], round[1])
    x[,2] <- base::round(x[,2], round[2])
  }
  character_columns <- which(apply(x, 2, function(column) {
    methods::is(tryCatch(as.numeric(column), warning = function(w) w), "warning")
  }))
  if(length(character_columns) > 0) {
    for(i in character_columns) {
      ## quote character columns
      x[,i] <- sprintf("\"%s\"", x[,1])
    }
  }
  paste(apply(x, 1, function(x) {
    paste(trimws(x), sep="", collapse=":")
  }), sep="", collapse=" ")
}
