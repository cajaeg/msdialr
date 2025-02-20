#' Export alignment results to MetaboAnalyst csv format
#'
#' @param x alignment results table
#' @param Sample vector of sample names
#' @param Class vector of class labels (needs to have >1 levels)
#' @param metabolite_names column name to use for metabolite names
#' @param file output file
#'
#' @returns NULL
#' @export
exportToMetaboAnalyst <- function(x, 
                                  Sample = attr(x, "msdial_sam")$name,
                                  Class = attr(x, "msdial_sam")$class,
                                  metabolite_names = "metabolite_name",
                                  file = "") {
  if(is.null(Sample))
    Sample <- formatNumericIDs(seq_along(getIntensityColumns(x)),
                               prefix = "Sample_")
  if(is.null(Class))
    Class <- formatNumericIDs(rep(1, length(getIntensityColumns(x))),
                              prefix = "Class_")
  if(length(unique(Class)) == 1)
    warning("MetaboAnalyst requires >1 factor level, 'Class' has only 1")
  if(metabolite_names %in% colnames(x))
    metabolite_names <- x[[ metabolite_names ]]
  else 
    metabolite_names <- formatNumericIDs(1:nrow(getIntensityMatrix(x, as.matrix = TRUE)),
                                         prefix = "Metabolite_")
  if(any(duplicated(metabolite_names)))
    metabolite_names <- make.unique(metabolite_names, "_")
  out <- data.frame(Sample = Sample,
                    Class = Class,
                    t(getIntensityMatrix(x, as.matrix = TRUE)))
  colnames(out)[3:ncol(out)] <- metabolite_names
  if(file != "") {
    message(
      'import this file on https://www.metaboanalyst.ca/MetaboAnalyst/upload/StatUploadView.xhtml',
      'using option "plain text file" -> "Peak intensities", "Samples in rows"'
    )
  }
  utils::write.csv(out, file = file)
}
