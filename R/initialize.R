#' Initial function
#'
#'
#' @param libname Path to package
#' @param pkgname Name of package
 .onLoad <- function(libname,pkgname) {
   libname=getwd()
   pkgname="SMRinR"
    if (Sys.info()["sysname"] == "Windows") {
      path=file.exists(paste0(libname,"/smr-1.3.1-win.exe"))
      if(!path){
        utils::download.file(url='https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-win-x86_64.zip',
                      destfile="smr.zip")
        files=c("smr-1.3.1-win-x86_64/libomp.dll",'smr-1.3.1-win-x86_64/libomp140.x86_64.dll',"smr-1.3.1-win-x86_64/smr-1.3.1-win.exe","smr-1.3.1-win-x86_64/zlib1.dll")
        utils::unzip("smr.zip",files = files,exdir=libname,junkpaths=T)}
    } else if (Sys.info()["sysname"] == "Linux") {
      path=file.exists(paste0(libname,"/smr"))
      if(!path){
        utils::download.file(url='https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip',
                      destfile="smr.zip")
        files=c("smr-1.3.1-linux-x86_64/smr")
        utils::unzip("smr.zip",files = files,exdir=libname,junkpaths=T)}
    }else {
      path=file.exists(paste0(libname,"/smr-1.3.1-macos-arm64"))
      if(!path){
        utils::download.file(url='https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip',
                      destfile="smr.zip")
        files=c("smr-1.3.1-macos-arm64")
        utils::unzip("smr.zip",files = files,exdir=libname,junkpaths=T)}
    }
  }
