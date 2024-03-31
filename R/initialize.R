
#' Initialize the R package of "SMRinR"
#'
#' @param file The path to the file of the R package initialization file
#' @examples
#' \donttest{
#' SMR_initialize('file path to initialization file')}
#' @export SMR_initialize
SMR_initialize=function(file){
  path <- system.file()
  unzip(file,overwrite = TRUE,exdir=gsub('base','SMRinR/tools',path))
}
