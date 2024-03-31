#' Get the path of SMR software
#'
#' @return the path of SMR software
#' @export
#'
smrfile=function(){
  if (Sys.info()["sysname"] == "Windows") {
    path <- system.file("tools/win", "smr.exe", package = "SMRinR", mustWork = TRUE)
  } else if (Sys.info()["sysname"] == "Linux") {
    path <- system.file("tools/linux", "smr", package = "SMRinR", mustWork = TRUE)
  } else {
    path <- system.file("tools/mac", "smr", package = "SMRinR", mustWork = TRUE)
  }
  return(path)
}

