#' Add chemical classification data from ClassyFire
#'
#' @param x MS-DIAL data.frame as returned by \code{\link{loadConsoleResults}()}
#'   or \code{\link{loadAlignmentResults}()}
#' @param in_column name of column containing InChIKeys
#' @param cache_file local cache file. Add
#'   \code{options(msdialr_classyfire_cache =
#'   "/path/to/classyfire_cache.sqlite")} (adjusted accordingly to your system)
#'   to your .Rprofile file. Use \code{Sys.getenv("R_USER")} to check where this
#'   file is located on your system.
#'
#' @return 'x' completed with the following columns: 'kingdom', 'superclass',
#'   'class', 'subclass', 'level5', 'level6'
#' @export
#' @importFrom rlang .data
#' @importFrom rlang :=
#'
#' @examples
#' \dontrun{
#' ik <- data.frame(inchikey = "MXWJVTOOROXGIU-UHFFFAOYSA-N")
#' res <- searchClassyFire(ik)
#' }
searchClassyFire <- function(x,
                             in_column = c("inchikey", "NISTres")[1],
                             cache_file = getOption("msdialr_classyfire_cache")
                             ) {
  if(!is.null(cache_file) && file.exists(cache_file)) {
    con <- classyfireR::open_cache(dbname = cache_file)
    on.exit(if(!is.null(con)) DBI::dbDisconnect(con))
  } else {
    con <- NULL
  }
  stopifnot(in_column %in% colnames(x))
  if (in_column == "NISTres") {
    ikeys <- do.call("rbind", x$NISTres)$inchikey
  } else {
    ikeys <- x[[ in_column ]]
  }
  ikeys <- ikeys[!is.na(ikeys)]
  ikeys <- ikeys[ikeys != ""]
  ikeys <- unique(ikeys)
  
  classy0 <- lapply(ikeys, classyfireR::get_classification, conn = con)
  classy1 <- lapply(classy0, function(x) {
    tryCatch(classyfireR::classification(x), error = function(e) NA)
  })
  names(classy1) <- sub(".*=", "", sapply(classy0, function(x) x@meta$inchikey))
  ClassyFireLevels <- c("kingdom", "superclass", "class", "subclass", "level5", "level6")
  classy <- plyr::ldply(classy1, .id = "inchikey", ) |>
    tibble::as_tibble() |> 
    dplyr::select(-"CHEMONT") |>
    dplyr::mutate("Level" := gsub(" ", "", .data$Level)) |>
    dplyr::filter(.data$Level %in% ClassyFireLevels) |>
    dplyr::mutate("Level" := factor(.data$Level, levels = ClassyFireLevels)) |>
    tidyr::pivot_wider(names_from = "Level", values_from = "Classification", names_expand = TRUE) |> 
    dplyr::mutate("inchikey" = as.character(levels(.data$inchikey))[.data$inchikey])
  x |> 
    dplyr::left_join(classy, by = "inchikey", suffix = c("", "$$")) |> 
    dplyr::select(-tidyselect::ends_with("$$"))
}
