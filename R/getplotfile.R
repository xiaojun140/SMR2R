
#' The function for get plot file
#'
#' @param ref  the parameter of ref only can be choosed one of the following: EUR, SAS, AFR, AMR, and EAS
#' @param gwas_path outcome data' path
#' @param qlt_path qlt data's path
#' @param probe The probe which you want to plot
#' @param hgfile specifies a gene range list(glist-hg18,glist-hg19 or glist-hg38).
#'
#' @return plot data
#' @export getplot
#'
getplot <- function(ref='EUR',
                   gwas_path,
                   qlt_path,
                   probe,
                   hgfile="glist-hg19"){
  ref=file.path(paste0(getwd(),"/","ref"))
  ref=paste0(ref,"/",list.files(ref,pattern = 'EUR')[1])
  ref=gsub("\\.[a-z]+$",'',ref)
  if(hgfile == "glist-hg18"){
    hgfile = paste0(find.package('SMR2R'),"/extdata/glist-hg18")}
  else if(hgfile == "glist-hg19"){
    hgfile = paste0(find.package('SMR2R'),"/extdata/glist-hg19")}
  else{
    hgfile = pahgfileste0(find.package('SMR2R'),"/extdata/glist-hg38")}

  command <- paste0(smrfile(),
                    " --bfile ",
                    ref,
                    " --gwas-summary ",
                    gwas_path,
                    " --beqtl-summary ",
                    qlt_path,
                    " --out ",
                    paste0(getwd(),"/myplot"),
                    " --plot",
                    " --probe ",
                    probe,
                    " --probe-wind 500 ",
                    "--gene-list ",
                    hgfile)
  system(command)
  SMRData = ReadSMRData(paste0(getwd(),"/plot/myplot.",probe,".txt"))
  return(SMRData)
}
