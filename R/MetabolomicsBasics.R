#' Wrapper for Metabolomics::RemoveFactorsByANOVA()
#'
#' @param x MS-DIAL alignment results table
#' @param sam sample table
#' @param fmod full ANOVA model
#' @param kmod reduced ANOVA model
#'
#' @returns an object of the same class as 'x'
#' @export
removeFactorsByANOVA <- function(x, 
                                 sam = NULL, 
                                 fmod = "class + injection_order",
                                 kmod = "class") {
  if(is.null(sam) && checkSam(x)) {
    sam <- attr(x, "msdial_sam")
    for(col_name in c("injection_order", "order")) {
      if(col_name %in% colnames(sam))
        sam[[ col_name ]] <- as.integer(sam[[ col_name ]])
    }
  }
  int_mat <- t(getIntensityMatrix(x, as.matrix = TRUE))
  int_mat_new <- 
    10^MetabolomicsBasics::RemoveFactorsByANOVA(
      log10(int_mat + 1),
      sam = sam,
      fmod = fmod,
      kmod = kmod,
      output = "y_norm"
    )
  x[, getIntensityColumns(x)] <- t(int_mat_new)
  x
}
