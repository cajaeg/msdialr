#' Add logP values to alignment results table
#'
#' @param x MS-DIAL alignment results table
#' @param in_column name of column containing InChIKey's. Can be "NISTres" in
#'   which case query is performed on a previous NIST search result.
#'
#' @return 'x' with two additional columns:
#' \itemize{
#' \item 'xlogp' XLogP value (see \code{rcdk::get.xlogp()})
#' \item 'alogp' ALogP value (see \code{rcdk::get.alogp()})
#' }
#' Existing columns are NOT replaced.
#' @export
#'
#' @examples
#' \dontrun{
#' x <- data.frame(smiles = c("C(C(CBr)Br)Cl", "CC1=C(C(=CC=C1)C)[N+](=O)[O-]"))
#' addLogP(x)
#' }
addLogP <- function(x, in_column = c("smiles", "NISTres")[1]) {
  stopifnot(in_column %in% colnames(x))
  if (in_column == "NISTres") {
    smiles <- do.call("rbind", x$NISTres)$smiles
  } else {
    smiles <- x[[ in_column ]]
  }
  smiles <- smiles[!is.na(smiles)]
  smiles <- smiles[smiles != ""]
  smiles <- unique(smiles)
  mol <- smi2mol(smiles)
  res <- data.frame(smiles = smiles, 
                    xlogp = NA_real_,
                    alogp = NA_real_)
  for (i in 1:length(mol)) {
    res$xlogp[i] <- rcdk::get.xlogp(mol[[i]])
    res$alogp[i] <- rcdk::get.alogp(mol[[i]])
  }
  out <- if (in_column == "NISTres") {
    x |>
      dplyr::mutate("NISTres" := purrr::map("NISTres", function(y) {
        dplyr::left_join(y, res, by = "smiles", suffix = c("", "$$")) |>
          dplyr::select(-tidyselect::ends_with("$$"))
      }))
  } else {
    x |>
      dplyr::left_join(res, by = "smiles", suffix = c("", "$$")) |>
      dplyr::select(-tidyselect::ends_with("$$"))
  }
  out |>
    relocateIntensityColumns() |>
    updateIntensityColumnIndex()
}


#' Convert SMILES to rcdk molecule
#'
#' @param smiles SMILES code (character vector)
#'
#' @return rcdk molecule
smi2mol <- function(smiles) {
  smiles_parser <- rcdk::get.smiles.parser()
  mol <- rcdk::parse.smiles(smiles, smiles.parser = smiles_parser)
  for(i in 1:length(mol)) {
    rcdk::set.atom.types(mol[[i]])
    rcdk::convert.implicit.to.explicit(mol[[i]])
    rcdk::do.isotopes(mol[[i]])
    rcdk::do.aromaticity(mol[[i]])
    mol[[i]] <- rcdk::generate.2d.coordinates(mol[[i]])
  }
  return(mol)
}


# addMolecularDescriptors <- function(x) {
#   # # descc <- rcdk::get.desc.categories()
#   # descn <- rcdk::get.desc.names()
#   # desc <- rcdk::eval.desc(mol, descn)
# }