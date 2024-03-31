
#' The function of running SMR in R
#'
#' @param ref the parameter of ref only can be choosed one of the following: EUR, SAS, AFR, AMR, and EAS
#' @param gwas_path outcome data' path
#' @param qlt_path qlt data's path
#' @param maf Numeric variable MAF (minor allele frequency)
#' @param out_file path of output file
#' @param thread the thread of work
#' @details
#' The demo code detailed in the test file
#'
#' @export runSMR
runSMR <- function(ref='EUR',#个体参考数据
                   gwas_path,#暴露数据
                   qlt_path,#qlt路径
                   maf=0.01,#maf值
                   out_file,#输出路径
                   thread = 4){
  ref=file.path(paste0(getwd(),"/","ref"))
  ref=paste0(ref,"/",list.files(ref,pattern = 'EUR')[1])
  ref=gsub("\\.[a-z]+$",'',ref)

  command <- paste0(smrfile(),
                    " --bfile ",
                    ref,
                    " --gwas-summary ",
                    gwas_path,
                    " --beqtl-summary ",
                    qlt_path,
                    " --maf ",
                    base::as.character(maf),
                    " --out ",
                    paste0(getwd(),'/',out_file)                    ,
                    paste0(" --thread-num ",base::as.character(thread)))
  system(command)
}

