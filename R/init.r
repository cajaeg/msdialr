.check_msdialr_options <- function() {
  msdialr_options <- list(
    msdialr_msdialconsole_path = c("MS-DIAL console", "/path/to/MSDIAL ver.4.9.221218 Windowsx64/MsdialConsoleApp.exe"),
    msdialr_msfinder_console = c("MS-FINDER console", "/path/to/MsfinderConsoleApp.exe"),
    msdialr_msconvert_path = c("MSConvert", "/path/to/msconvert.exe"),
    msdialr_mspepsearch_path = c("MSPepSearch", "/path/to/MSPepSearch.exe"),
    msdialr_classyfire_cache = c("ClassyFire", "/path/to/classyfire_cache.rda")
    # msdialr_pubchem_sqlite = "/path/to/pubchem_sqlite.rda",
    # "msdialr_localNISTrda" = "/path/to/NIST.rda",
  )
  for (i in 1:length(msdialr_options)) {
    if (is.null(getOption(names(msdialr_options)[i]))) {
      packageStartupMessage(sprintf(
        "Hint: to use %s, add 'options(%s = \"%s\")' to your .Rprofile",
        msdialr_options[[i]][1], names(msdialr_options)[i], msdialr_options[[i]][2]
      ))
    } else {
      if(!file.exists(getOption(names(msdialr_options)[i]))) 
        packageStartupMessage(sprintf(
          "Warning: option '%s' does not point to valid file",
          names(msdialr_options)[i]
        ))
    }
  }
}