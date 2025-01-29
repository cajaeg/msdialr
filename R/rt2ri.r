#' Calculate RT from RI based on linear interpolation
#'
#' @param ri numeric()
#' @param ref alkane data.frame as returned by \code{\link{loadAlkanes}()}
#'   (Num|RT|RI)
#' @return numeric()
#' @export
#'
#' @examples
#' fp <- system.file("extdata/Alkane.txt", package = "msdialr")
#' alk <- loadAlkanes(fp)
#' head(alk)
#' ri2rt(ri = c(1000,1200,1400), ref = alk)
ri2rt <- function(ri, ref) {
  stopifnot(all(diff(ref$RT) > 0))
  do_it <- function(ri, ref) {
    if (is.finite(ri)) {
      if (ri < min(ref$RI)) {
        b <- c(1, 2)
      } else if (ri >= max(ref$RI)) {
        b <- nrow(ref) - c(1, 0)
      } else {
        b <- c(max(which(ref$RI <= ri)), min(which(ref$RI > ri)))
      }
      fit <- stats::lm(RT ~ RI, data = ref[b, ])
      as.numeric(stats::predict(fit, newdata = data.frame(RI = ri)))
    } else {
      NA_real_
    }
  }
  sapply(ri, do_it, ref = ref)
}


#' Calculate RI from RT based on linear interpolation
#'
#' @param rt numeric()
#' @param ref alkane data.frame as returned by \code{\link{loadAlkanes}()}
#'   (Num|RT|RI)
#' @return numeric()
#' @export
#'
#' @examples
#' fp <- system.file("extdata/Alkane.txt", package = "msdialr")
#' alk <- loadAlkanes(fp)
#' head(alk)
#' rt2ri(rt = c(5,10,15), ref = alk)
rt2ri <- function(rt, ref) {
  stopifnot(all(diff(ref$RT) > 0))
  do_it <- function(rt, ref) {
    if (is.finite(rt)) {
      if (rt < min(ref$RT)) {
        b <- c(1, 2)
      } else if (rt >= max(ref$RT)) {
        b <- nrow(ref) - c(1, 0)
      } else {
        b <- c(max(which(ref$RT <= rt)), min(which(ref$RT > rt)))
      }
      fit <- stats::lm(RI ~ RT, data = ref[b, ])
      as.numeric(stats::predict(fit, newdata = data.frame(RT = rt)))
    } else {
      NA_real_
    }
  }
  sapply(rt, do_it, ref = ref)
}

