#' Get the path of SMR software
#'
#' @return the path of SMR software
#' @export
#'
smrfile=function(){
  if (Sys.info()["sysname"] == "Windows") {
    path=paste0(base::find.package("SMRinR"),"/inst/extdata","/smr-1.3.1-win.exe")
  } else if (Sys.info()["sysname"] == "Linux") {
    path=paste0(base::find.package("SMRinR"),"/inst/extdata","/smr")
  } else {
    path=paste0(base::find.package("SMRinR"),"/inst/extdata","/smr-1.3.1-macos-arm64")
  }
  return(path)
}

