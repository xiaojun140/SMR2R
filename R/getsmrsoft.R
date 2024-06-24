#' Get the path of SMR software
#'
#' @return the path of SMR software
#' @export
#'
smrfile=function(){
  initialize()
  if (Sys.info()["sysname"] == "Windows") {
    path=paste0(getwd(),"/smr-1.3.1-win.exe")
  } else if (Sys.info()["sysname"] == "Linux") {
    path=paste0(getwd(),"/inst/extdata","/smr")
  } else {
    path=paste0(getwd(),"/inst/extdata","/smr-1.3.1-macos-arm64")
  }
  return(path)
}

