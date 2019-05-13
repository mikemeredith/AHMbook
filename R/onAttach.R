

.onAttach <- function(libname, pkgname) {
  version <- try(packageVersion('AHMbook'), silent=TRUE)
  if(!inherits(version, "try-error"))
    packageStartupMessage("This is AHMbook ", version,
      ". For overview type ?AHMbook; for changes do news(p='AHMbook').")
}
