#' Initial function
#'
#'
#' @param libname Path to package
#' @param pkgname Name of package
 .onLoad <- function(libname, pkgname) {
   libname=base::find.package("SMRinR")
   pkgname="SMRinR"
    if (Sys.info()["sysname"] == "Windows") {
      path=file.exists(paste0(libname,"/inst/extdata","/smr-1.3.1-win.exe"))
      if(!path){
        utils::download.file(url='https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-win-x86_64.zip',
                      destfile="smr.zip")
        files=c("smr-1.3.1-win-x86_64/libomp.dll",'smr-1.3.1-win-x86_64/libomp140.x86_64.dll',"smr-1.3.1-win-x86_64/smr-1.3.1-win.exe","smr-1.3.1-win-x86_64/zlib1.dll")
        utils::unzip("smr.zip",files = files,exdir=paste0(find.package('SMRinR'),"/inst/extdata"),junkpaths=T)}
    } else if (Sys.info()["sysname"] == "Linux") {
      path=file.exists(paste0(libname,"/inst/extdata","/smr"))
      if(!path){
        utils::download.file(url='https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip',
                      destfile="smr.zip")
        files=c("smr")
        utils::unzip("smr.zip",files = files,exdir=paste0(find.package('SMRinR'),"/inst/extdata"),junkpaths=T)}
    }else {
      path=file.exists(paste0(libname,"/inst/extdata","/smr-1.3.1-macos-arm64"))
      if(!path){
        utils::download.file(url='https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip',
                      destfile="smr.zip")
        files=c("smr-1.3.1-macos-arm64")
        utils::unzip("smr.zip",files = files,exdir=paste0(find.package('SMRinR'),"/inst/extdata"),junkpaths=T)}
    }

  }
